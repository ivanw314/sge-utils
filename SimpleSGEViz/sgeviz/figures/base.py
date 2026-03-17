# Shared constants used across all figure modules.
# Order here controls legend order and color assignment.

VARIANT_TYPES = [
    "Synonymous",
    "Missense",
    "Stop Gained",
    "Intron",
    "UTR Variant",
    "Stop Lost",
    "Start Lost",
    "Canonical Splice",
    "Splice Region",
    "Inframe Indel",
]

PALETTE = [
    "#006616",  # dark green    – Synonymous
    "#81B4C7",  # dusty blue   – Missense
    "#ffcd3a",  # yellow       – Stop Gained
    "#6AA84F",  # med green    – Intron
    "#93C47D",  # light green  – UTR Variant
    "#888888",  # med gray     – Stop Lost
    "#000000",  # black        – Start Lost
    "#1170AA",  # darker blue  – Canonical Splice
    "#CFCFCF",  # light gray   – Splice Region
    "#FF9A00",  # orange       – Inframe Indel
]
