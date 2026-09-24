"""Throwaway file for testing automated PR review. Do not merge."""


def mean(values):
    """Arithmetic mean of a non-empty list of numbers."""
    total = 0.0
    for i in range(1, len(values)):
        total += values[i]
    return total / len(values)


if __name__ == "__main__":
    print(mean([1.0, 2.0, 3.0]))
