from __future__ import annotations

from collections import Counter
import math
import random
from collections.abc import Sequence


ALLOWED_ALPHABET = ("A", "C", "G", "T")
ALLOWED_PAIRS = (
    ("A", "T"),
    ("C", "G"),
    ("G", "T"),
    ("T", "A"),
    ("G", "C"),
    ("T", "G"),
)

_PATTERN_TRANSLATION = str.maketrans({"A": "T", "C": "G", "G": "Y", "T": "R"})
_COMPATIBLE = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
}

Coord = tuple[int, int]


def _normalize_sequence(seq: str | Sequence[str]) -> str:
    if isinstance(seq, str):
        return seq.upper()
    return "".join(str(base) for base in seq).upper()


def _has_ambiguous_bases(seq: str) -> bool:
    return any(base not in ALLOWED_ALPHABET for base in seq)


def _pattern_sequence(seq: str) -> str:
    return seq.translate(_PATTERN_TRANSLATION)


def _matches_at(pattern: str, text: str, start: int) -> bool:
    for offset, pattern_base in enumerate(pattern):
        if text[start + offset] not in _COMPATIBLE.get(pattern_base, ()):
            return False
    return True


def _find_all(pattern: str, text: str) -> list[Coord]:
    pattern_len = len(pattern)
    if pattern_len == 0 or pattern_len > len(text):
        return []
    return [
        (idx + 1, idx + pattern_len)
        for idx in range(len(text) - pattern_len + 1)
        if _matches_at(pattern, text, idx)
    ]


def _find_first(pattern: str, text: str) -> Coord | None:
    matches = _find_all(pattern, text)
    return matches[0] if matches else None


def _find_last(pattern: str, text: str) -> Coord | None:
    matches = _find_all(pattern, text)
    return matches[-1] if matches else None


def _occurs_in(pattern: str, text: str) -> bool:
    return _find_first(pattern, text) is not None


def _shift(coord: Coord, offset: int) -> Coord:
    return coord[0] + offset, coord[1] + offset


def _coord_len(coord: Coord) -> int:
    return coord[1] - coord[0] + 1


def _coerce_coord_range(seq_range: Coord | range | Sequence[int]) -> Coord:
    if isinstance(seq_range, range):
        if seq_range.step != 1:
            raise ValueError("seqA_range must have a step of 1")
        if len(seq_range) == 0:
            raise ValueError("seqA_range cannot be empty")
        return seq_range.start, seq_range.stop - 1
    if len(seq_range) != 2:
        raise ValueError("seqA_range must be a (start, end) pair or a Python range")
    start, end = int(seq_range[0]), int(seq_range[1])
    if start > end:
        raise ValueError("seqA_range start must be <= end")
    return start, end


def comp_segm_len_to_ds_force(
    csl: int,
    alpha: float,
    seqlen: int,
    seqlen_forced: int | None = None,
    c0: float = -2.2,
) -> float:
    """Convert the longest complementary segment length to DS force."""
    if csl < 2:
        return math.nan
    if seqlen_forced is None:
        seqlen_forced = seqlen
    return math.log(1 / alpha) - (math.log(seqlen) + math.log(seqlen_forced)) / (csl - c0)


def compute_alpha(seq: str | Sequence[str]) -> float:
    """Compute the probability that two randomly chosen bases form an allowed pair."""
    seq = _normalize_sequence(seq)
    seqlen = len(seq)
    nt_usage = Counter(seq)
    if sum(nt_usage.get(nt, 0) for nt in ALLOWED_ALPHABET) != seqlen:
        return 0.0
    return sum(nt_usage[a] * nt_usage[b] / seqlen**2 for a, b in ALLOWED_PAIRS)


def find_longest_complementary_segment(
    seq: str | Sequence[str],
    seq_a_range: Coord | range | Sequence[int] | None = None,
    return_coords: bool = False,
    *,
    rand_choice: bool = False,
    seed: int | None = None,
) -> int | tuple[Coord, Coord]:
    """Find the longest compatible double-stranded segment.

    Returned coordinates are 1-based inclusive pairs.
    """
    seq = _normalize_sequence(seq)
    if seq_a_range is not None:
        return _find_longest_complementary_segment_constrained(seq, seq_a_range, return_coords)

    pattern_seq = _pattern_sequence(seq)
    seq_len = len(seq)
    max_len = 2
    start_possibilities = [1]

    for i in range(1, seq_len + 1):
        for k in range(i + max_len - 1, seq_len + 1):
            pattern = pattern_seq[i - 1 : k][::-1]
            if not _occurs_in(pattern, seq[k:]):
                if (k - i) > max_len:
                    max_len = k - i
                    start_possibilities = [i]
                elif (k - i) == max_len:
                    start_possibilities.append(i)
                break

    if not return_coords:
        return max_len

    rng = random.Random(seed)
    start_pos = rng.choice(start_possibilities) if rand_choice else start_possibilities[-1]
    pattern = pattern_seq[start_pos - 1 : start_pos + max_len - 1][::-1]
    downstream_possibilities = _find_all(pattern, seq[start_pos + max_len - 1 :])
    if not downstream_possibilities:
        return (seq_len, seq_len), (seq_len, seq_len)

    downstream = rng.choice(downstream_possibilities) if rand_choice else downstream_possibilities[0]
    return (start_pos, start_pos + max_len - 1), _shift(downstream, start_pos + max_len - 1)


def _find_longest_complementary_segment_constrained(
    seq: str,
    seq_a_range: Coord | range | Sequence[int],
    return_coords: bool = False,
) -> int | tuple[Coord, Coord]:
    pattern_seq = _pattern_sequence(seq)
    seq_len = len(seq)
    range_start, range_end = _coerce_coord_range(seq_a_range)
    if range_start < 1 or range_end > seq_len:
        raise ValueError("seqA_range must be within the sequence coordinates")

    max_len = 2
    start_pos = range_start
    for i in range(range_start, range_end + 1):
        for k in range(i + max_len - 1, range_end + 1):
            pattern = pattern_seq[i - 1 : k][::-1]
            has_pre_match = _occurs_in(pattern, seq[: i - 1])
            has_post_match = _occurs_in(pattern, seq[k:])
            if not has_pre_match and not has_post_match:
                if (k - i) > max_len:
                    max_len = k - i
                    start_pos = i
                elif (k - i) == max_len:
                    start_pos = i
                break
            if k == range_end:
                if (k - i + 1) > max_len:
                    max_len = k - i + 1
                    start_pos = i
                elif (k - i + 1) == max_len:
                    start_pos = i

    if not return_coords:
        return max_len

    pattern = pattern_seq[start_pos - 1 : start_pos + max_len - 1][::-1]
    pre_match = _find_last(pattern, seq[: start_pos - max_len])
    post_match = _find_first(pattern, seq[start_pos + max_len - 1 :])
    if pre_match is None and post_match is None:
        return (seq_len, seq_len), (seq_len, seq_len)
    if post_match is not None:
        return (start_pos, start_pos + max_len - 1), _shift(post_match, start_pos + max_len - 1)
    return pre_match, (start_pos, start_pos + max_len - 1)


def find_longest_complementary_segment_last(seq: str | Sequence[str]) -> tuple[Coord, Coord]:
    """Find the longest complementary segment ending at the end of the sequence."""
    seq = _normalize_sequence(seq)
    pattern_seq = _pattern_sequence(seq)
    seq_len = len(seq)
    max_len = 2
    start_pos = seq_len - 2

    for k in range(max_len, seq_len + 1):
        pattern = pattern_seq[seq_len - k - 1 : seq_len][::-1]
        if not _occurs_in(pattern, seq[: seq_len - k - 1]):
            if k + 1 >= max_len:
                max_len = k + 1
                start_pos = seq_len - k
            break

    pattern = pattern_seq[start_pos:seq_len][::-1]
    first_match = _find_first(pattern, seq[:start_pos])
    if first_match is None:
        return (seq_len, seq_len), (seq_len, seq_len)
    return first_match, (start_pos + 1, seq_len)


def find_longest_complementary_segment_sliding(
    seq: str | Sequence[str],
    sliding_window_length: int | None = None,
) -> tuple[list[int], list[float], list[Coord], list[Coord]]:
    """Run the complementary-segment search on stride-1 sliding windows."""
    seq = _normalize_sequence(seq)
    if sliding_window_length is None:
        sliding_window_length = min(len(seq), 3000)
    if sliding_window_length < 1:
        raise ValueError("sliding_window_length must be positive")

    csl_values: list[int] = []
    alphas: list[float] = []
    range1s: list[Coord] = []
    range2s: list[Coord] = []
    global_pos = 1

    while global_pos <= len(seq) - sliding_window_length + 1:
        sliding_seq = seq[global_pos - 1 : global_pos + sliding_window_length - 1]
        alpha = compute_alpha(sliding_seq)
        if alpha == 0:
            csl_values.append(0)
            alphas.append(0.0)
            range1s.append((1, 1))
            range2s.append((1, 1))
            global_pos += 1
            continue

        alphas.append(alpha)
        if global_pos == 1:
            range1, range2 = find_longest_complementary_segment(sliding_seq, return_coords=True)
            range1s.append(_shift(range1, global_pos - 1))
            range2s.append(_shift(range2, global_pos - 1))
        else:
            previous_range1 = range1s[-1]
            if previous_range1[0] < global_pos:
                range1, range2 = find_longest_complementary_segment(sliding_seq, return_coords=True)
                range1s.append(_shift(range1, global_pos - 1))
                range2s.append(_shift(range2, global_pos - 1))
            else:
                previous_length = csl_values[-1]
                last_range1, last_range2 = find_longest_complementary_segment_last(sliding_seq)
                if _coord_len(last_range1) > previous_length:
                    range1s.append(_shift(last_range1, global_pos - 1))
                    range2s.append(_shift(last_range2, global_pos - 1))
                else:
                    range1s.append(range1s[-1])
                    range2s.append(range2s[-1])

        csl_values.append(_coord_len(range1s[-1]))
        global_pos += 1

    return csl_values, alphas, range1s, range2s


def compute_ds_force(
    seq: str | Sequence[str],
    *,
    return_lcs_positions: bool = False,
    sliding_window_length: int | None = None,
    seq_a_range: Coord | range | Sequence[int] | None = None,
) -> float | None | list[float] | tuple[float, Coord, Coord] | tuple[list[float], list[Coord], list[Coord]]:
    """Compute DS force for a sequence.

    Ambiguous bases return ``None`` to mirror Julia's ``missing``.
    Coordinates are returned as 1-based inclusive ``(start, end)`` tuples.
    """
    seq = _normalize_sequence(seq)
    if _has_ambiguous_bases(seq):
        return None

    if sliding_window_length is not None:
        csl_values, alphas, range1s, range2s = find_longest_complementary_segment_sliding(
            seq, sliding_window_length
        )
        ds_forces = [
            comp_segm_len_to_ds_force(csl, alpha, sliding_window_length)
            for csl, alpha in zip(csl_values, alphas, strict=True)
        ]
        if return_lcs_positions:
            return ds_forces, range1s, range2s
        return ds_forces

    alpha = compute_alpha(seq)
    if return_lcs_positions:
        range_a, range_b = find_longest_complementary_segment(
            seq, seq_a_range=seq_a_range, return_coords=True
        )
        csl = _coord_len(range_a)
        if seq_a_range is None:
            return comp_segm_len_to_ds_force(csl, alpha, len(seq)), range_a, range_b
        forced_start, forced_end = _coerce_coord_range(seq_a_range)
        return (
            comp_segm_len_to_ds_force(csl, alpha, len(seq), forced_end - forced_start + 1),
            range_a,
            range_b,
        )

    csl = find_longest_complementary_segment(seq, seq_a_range=seq_a_range)
    if seq_a_range is None:
        return comp_segm_len_to_ds_force(csl, alpha, len(seq))
    forced_start, forced_end = _coerce_coord_range(seq_a_range)
    return comp_segm_len_to_ds_force(csl, alpha, len(seq), forced_end - forced_start + 1)


CompSegmLenToDSForce = comp_segm_len_to_ds_force
ComputeAlpha = compute_alpha
ComputeDSForce = compute_ds_force
FindLongestComplementarySegment = find_longest_complementary_segment
FindLongestComplementarySegmentLast = find_longest_complementary_segment_last
FindLongestComplementarySegmentSliding = find_longest_complementary_segment_sliding
