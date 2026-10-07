# %%
import json
import pandas as pd

# %%
with open("./components.json", "r") as f:
    components = json.load(f)

with open("./formulation_vals.json") as f:
    formulation_vals = json.load(f)
# %%

# Constants
MRNA_LENGTH = 4000
TOTAL_MRNA = 5000
NP_RATIO = 6
MASTER_MIX = 1.6
MRNA_CONC = 1000
VOL_RATIO = 3
MRNA_CONC_IN_LNP = 112.5


# Leave the formulas unsimplified for Mae to track things
# %%


class LNPCalculator:

    def __init__(self, total_mrna, mrna_conc_in_lnp, vol_ratio):
        self.vol_ratio = vol_ratio
        self.total_mrna = total_mrna
        self.mrna_conc_in_lnp = mrna_conc_in_lnp
        self.total_vol_mrna_lnp = self._total_vol_mrna_lnp(
            self.total_mrna, self.mrna_conc_in_lnp
        )
        self.total_vol_lipids = self._total_vol_lipids(
            self.vol_ratio, self.total_vol_mrna_lnp
        )
        self.total_vol_mrna = self._total_mrna_vol(
            self.vol_ratio, self.total_vol_mrna_lnp
        )

        return self

    def p_ratio(total_mrna: int, mrna_length: int) -> float:
        ratio = ((total_mrna * 1e-9) / (mrna_length * 320.5 + 159)) * 1e12 * mrna_length
        return ratio

    def vol_a(
        mol_weight_a: float,
        np_ratio: float,
        p_ratio: float,
        master_mix: float,
        stock_conc_a: float,
    ) -> float:
        vol_for_a = (
            (mol_weight_a * 1e6 * np_ratio * p_ratio * 1e-12)
            / 10
            * master_mix
            * 10
            / stock_conc_a
        )
        return vol_for_a

    def vol_b(
        mol_ratio_b,
        mol_ratio_a,
        p_ratio,
        np_ratio,
        mol_weight_b,
        master_mix,
        stock_conc_b,
    ) -> float:

        num = (
            mol_ratio_b / mol_ratio_a * p_ratio * np_ratio * 1e-12 * mol_weight_b * 1e6
        )

        vol_for_b = (num / 10) * master_mix * (10 / stock_conc_b)

        return vol_for_b

    def _total_vol_mrna_lnp(total_mrna, mrna_conc_in_lnp):
        return total_mrna / mrna_conc_in_lnp

    def _total_vol_lipids(vol_ratio, total_vol_mrna_lnp):
        return total_vol_mrna_lnp / (1 + vol_ratio)

    def _total_mrna_vol(vol_ratio, total_vol_mrna_lnp):
        return total_vol_mrna_lnp * vol_ratio / (vol_ratio + 1)

    def citrate_buffer_mul(total_mrna_vol):
        return total_mrna_vol / 10

    def mrna_vol_mul(total_mrna, mrna_conc):
        return total_mrna / mrna_conc

    def water_remainder(total_mrna_vol, total_mrna, mrna_conc):
        cit_buf = self.citrate_buffer_ml(total_mrna)
        mrna_vol = self.mrna_vol_ml(total_mrna, mrna_conc)

        return total_mrna_vol - cit_buf - mrna_vol

    def lipid_to_mrna(total_vol_mrna_lnp, vol_ratio):
        return total_vol_mrna_lnp / (vol_ratio + 1)

    def ethanol_mul(total_vol_lipids, master_mix, vol_components: dict):

        sum_vol_component = 0
        for key, val in vol_components:
            vol_comp = self.vol_b()

        return
