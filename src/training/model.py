"""Weaver --network-config: three-stream ParT (pfcand + v0njet + sv), full
architecture sized for a real GPU (e.g. NVIDIA A100) rather than the M1's
8-core GPU -- pair_embed_dims stays at [64, 64, 64] here (vs. the M1
variant's [32, 32, 32]). This is the 5-class B/C/S/D/U version.

weaver's built-in `ParticleTransformerTagger` (in weaver.nn.model.ParticleTransformer)
only combines TWO object streams (pf + one other, normally "sv"). To add a
THIRD stream (v0njet) alongside pf and sv, this file defines a small custom
module that mirrors ParticleTransformerTagger's own internal pattern, just
extended from two Embed+Trimmer pairs to three, all feeding into one shared
`ParticleTransformer` core -- built entirely from that module's public
classes (Embed, SequenceTrimmer, ParticleTransformer), no changes to the
library itself.

As with the two-stream version: weaver calls `model(*inputs)` with `inputs`
built from `data_config.input_names` in declaration order, so this wrapper's
forward() signature must positionally match data_config_ALEPH_ParT_full.yaml's
`inputs:` block exactly:
    pf_features, pf_vectors, pf_mask, v0_features, v0_vectors, v0_mask,
    sv_features, sv_vectors, sv_mask
Reorder one, reorder the other, or a tensor silently lands in the wrong slot.

No device-specific code needed for CUDA: device selection happens entirely
via weaver's --gpus flag, and this file makes no MPS/CUDA-specific
assumptions -- it runs correctly on either backend as-is.
"""

import torch
from torch import nn

from weaver.nn.model.ParticleTransformer import Embed, ParticleTransformer, SequenceTrimmer


class ParticleTransformerThreeStreamTagger(nn.Module):
    def __init__(
        self,
        pf_input_dim,
        v0_input_dim,
        sv_input_dim,
        num_classes=None,
        pair_input_type="pp",
        pair_input_dim=4,
        pair_extra_dim=0,
        remove_self_pair=False,
        use_pre_activation_pair=True,
        embed_dims=(128, 512, 128),
        pair_embed_dims=(64, 64, 64),
        num_heads=8,
        num_layers=8,
        num_cls_layers=2,
        block_params=None,
        cls_block_params=None,
        fc_params=(),
        activation="gelu",
        version=1,
        weight_init="moco",
        fix_init=True,
        trim=True,
        for_inference=False,
        for_segmentation=False,
        use_amp=False,
        **kwargs,
    ):
        super().__init__()
        self.use_amp = use_amp

        # one trimmer + one embed per stream, same pattern ParticleTransformerTagger
        # uses for its two streams (see weaver/nn/model/ParticleTransformer.py)
        self.pf_trimmer = SequenceTrimmer(enabled=trim and not for_inference)
        self.v0_trimmer = SequenceTrimmer(enabled=trim and not for_inference)
        self.sv_trimmer = SequenceTrimmer(enabled=trim and not for_inference)

        self.pf_embed = Embed(pf_input_dim, embed_dims, activation=activation)
        self.v0_embed = Embed(v0_input_dim, embed_dims, activation=activation)
        self.sv_embed = Embed(sv_input_dim, embed_dims, activation=activation)

        # single shared core: takes the concatenated (pf + v0 + sv) sequence,
        # so pairwise-interaction features are computed across ALL THREE
        # streams together, not per-stream.
        self.part = ParticleTransformer(
            input_dim=embed_dims[-1],
            num_classes=num_classes,
            pair_input_type=pair_input_type,
            pair_input_dim=pair_input_dim,
            pair_extra_dim=pair_extra_dim,
            remove_self_pair=remove_self_pair,
            use_pre_activation_pair=use_pre_activation_pair,
            embed_dims=[],  # already embedded per-stream above
            pair_embed_dims=pair_embed_dims,
            num_heads=num_heads,
            num_layers=num_layers,
            num_cls_layers=num_cls_layers,
            block_params=block_params,
            cls_block_params=cls_block_params,
            fc_params=fc_params,
            activation=activation,
            version=version,
            weight_init=weight_init,
            fix_init=fix_init,
            trim=False,  # trimming is handled per-stream above, before concatenation
            for_inference=for_inference,
            for_segmentation=for_segmentation,
            use_amp=use_amp,
        )

    @torch.jit.ignore
    def no_weight_decay(self):
        return {"part.cls_token"}

    def forward(self, pf_x, pf_v, pf_mask, v0_x, v0_v, v0_mask, sv_x, sv_v, sv_mask):
        # x: (N, C, P)   v: (N, 4, P) [px,py,pz,energy]   mask: (N, 1, P)
        with torch.no_grad():
            pf_x, pf_v, pf_mask, _ = self.pf_trimmer(pf_x, pf_v, pf_mask)
            v0_x, v0_v, v0_mask, _ = self.v0_trimmer(v0_x, v0_v, v0_mask)
            sv_x, sv_v, sv_mask, _ = self.sv_trimmer(sv_x, sv_v, sv_mask)
            v = torch.cat([pf_v, v0_v, sv_v], dim=2)
            mask = torch.cat([pf_mask, v0_mask, sv_mask], dim=2)

        pf_x = self.pf_embed(pf_x)  # -> (batch, seq_len, embed_dim)
        v0_x = self.v0_embed(v0_x)
        sv_x = self.sv_embed(sv_x)
        x = torch.cat([pf_x, v0_x, sv_x], dim=1)

        return self.part(x, v, mask)


def get_model(data_config, **kwargs):
    pf_input_dim = len(data_config.input_dicts["pf_features"])
    v0_input_dim = len(data_config.input_dicts["v0_features"])
    sv_input_dim = len(data_config.input_dicts["sv_features"])
    num_classes = len(data_config.label_value)

    cfg = dict(
        pf_input_dim=pf_input_dim,
        v0_input_dim=v0_input_dim,
        sv_input_dim=sv_input_dim,
        num_classes=num_classes,
        pair_input_type="pp",
        pair_input_dim=4,  # deltaR, kT, z, m^2
        remove_self_pair=False,
        use_pre_activation_pair=True,
        embed_dims=[128, 512, 128],
        pair_embed_dims=[64, 64, 64],
        num_heads=8,
        num_layers=8,
        num_cls_layers=2,
        block_params=None,
        cls_block_params={"dropout": 0, "attn_dropout": 0, "activation_dropout": 0},
        fc_params=[(128, 0.1)],
        activation="gelu",
        version=1,  # bump to 2/3 for the SwiGLU/RMSNorm variants discussed earlier
        weight_init="moco",
        fix_init=True,
        trim=True,
        for_inference=False,
        for_segmentation=False,
    )
    cfg.update(kwargs)  # lets --network-option override anything above from the CLI

    model = ParticleTransformerThreeStreamTagger(**cfg)
    #model = model.double()
    model_info = {
        "input_names": list(data_config.input_names),
        "input_shapes": {k: ((1,) + s[1:]) for k, s in data_config.input_shapes.items()},
        "output_names": ["softmax"],
        "dynamic_axes": {
            **{k: {0: "N", 2: "n_" + k.split("_")[0]} for k in data_config.input_names},
            "softmax": {0: "N"},
        },
    }

    return model, model_info


def get_loss(data_config, **kwargs):
    return torch.nn.CrossEntropyLoss()
