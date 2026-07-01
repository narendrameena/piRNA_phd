"""
Canonical strain order for EVERY figure in this project.

Source: thesis Figure 4.4 (RNA-seq PCA) strain legend; strains ordered by median
P20.5 PC1 position from the RNA-seq PCA (thesis p.136). Same order used in the
piRNA chapter (Fig 5.7/5.9). This is NOT the Snakemake/project order and NOT
alphabetical. Import this everywhere instead of redefining order per script.

  from strain_order import STRAIN_ORDER, WILD, TIMEPOINT_ORDER, strain_rank
"""

STRAIN_ORDER = [
    'C57BL_6', 'C57BL_6NJ', 'BALB_cJ', 'A_J', 'FVB_NJ', 'C3H_HeJ', 'LP_J',
    '129S1_SvImJ', 'DBA_2J', 'AKR_J', 'CBA_J', 'NZO_HlLtJ', 'NOD_ShiLtJ',
    'WSB_EiJ', 'CAST_EiJ', 'PWK_PhJ', 'SPRET_EiJ',
]

# C57BL_6 (pos 1) = reference C57BL/6J (T2T), distinct from C57BL_6NJ.

WILD = {'CAST_EiJ', 'PWK_PhJ', 'SPRET_EiJ', 'WSB_EiJ'}

TIMEPOINT_ORDER = ['E16.5', 'P12.5', 'P20.5']  # developmental, never alphabetical

_RANK = {s: i for i, s in enumerate(STRAIN_ORDER)}

def strain_rank(s):
    """Sort key: position in STRAIN_ORDER; unknown strains sort to the end."""
    return _RANK.get(s, len(STRAIN_ORDER))
