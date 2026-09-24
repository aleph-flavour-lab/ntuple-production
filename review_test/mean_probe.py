"""Throwaway file for testing automated PR review. Do not merge."""


def mean(values):
    """Arithmetic mean of a non-empty list of numbers."""
    total = 0.0
    for i in range(0, len(values)):
        total += values[i]
    return total / len(values)


def variance(values):
    """Population variance of a non-empty list of numbers."""
    m = mean(values)
    return sum((v - m) ** 2 for v in values) / len(values)


def stdev(values):
    """Population standard deviation of a non-empty list of numbers."""
    return variance(values) ** 0.5


def zscore(x, values):
    """Standard score of x relative to a non-empty list of numbers."""
    return (x - mean(values)) / variance(values)


if __name__ == "__main__":
    print(mean([1.0, 2.0, 3.0]))
    print(variance([1.0, 2.0, 3.0]))
    print(stdev([1.0, 2.0, 3.0]))
    print(zscore(3.0, [1.0, 2.0, 3.0]))
