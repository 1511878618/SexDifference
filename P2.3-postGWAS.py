#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@Description:       :
@Date     :2025/01/13 16:27:24
@Author      :Tingfeng Xu
@version      :1.0
"""
import argparse
import json
import textwrap
from pathlib import Path

import gwaslab as gl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gwaslab.g_Log import Log
from scipy import stats
from scipy.stats import norm


def cal_p_from_z(z):
    return 2 * norm.sf(np.abs(z))


def calculate_lambda_from_z(z_scores):
    chisq = z_scores**2
    median_chisq = np.median(chisq)
    lambda_gc = median_chisq / 0.4549364
    return lambda_gc


def calculate_lambda_from_p(p_values):
    z_scores = np.abs(stats.norm.ppf(p_values / 2))
    return calculate_lambda_from_z(z_scores)


def save_fig(
    fig=None,
    path=None,
    bbox_inches="tight",
    dpi=400,
    tiff=False,
    pdf=False,
    svg=False,
    tiff_compress=False,
    **kwargs,
):
    if path is None:
        path = "temp"
    if fig is None:
        fig = plt.gcf()
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(f"{path}.png", dpi=dpi, bbox_inches=bbox_inches, **kwargs)
    if pdf:
        fig.savefig(f"{path}.pdf", dpi=dpi, bbox_inches=bbox_inches, **kwargs)
    if svg:
        fig.savefig(f"{path}.svg", dpi=dpi, bbox_inches=bbox_inches, **kwargs)
    if tiff:
        fig.savefig(f"{path}.tiff", dpi=dpi, bbox_inches=bbox_inches, **kwargs)
    # plt.close(fig)


def het_test(BETA_F, SE_F, BETA_M, SE_M, P_F=None, P_M=None, r=None):
    """
    Follow the formula in the paper of (Bernabeu et al., Nat Genet, 2021)

    BETA_F: np.array or pd.Series,
    SE_F: np.array or pd.Series,
    BETA_M: np.array or pd.Series,
    SE_M: np.array or pd.Series,
    r: float, default None, the spearman correlation between the Female and Male of p-values (zscore is same as p-value)
    """

    if r is None:
        # r, _ = spearmanr(BETA_F, BETA_M)
        if P_F is None:
            P_F = 2 * norm.sf(np.abs(BETA_F / SE_F))
        if P_M is None:
            P_M = 2 * norm.sf(np.abs(BETA_M / SE_M))

        r, p_spearmanr = spearmanr(P_F, P_M)

        print(
            f"The spearman correlation of Z-scores is {r:.6f}, p-value is {p_spearmanr:.4e}"
        )
    # calculate the two-tailed t-test
    t = (BETA_M - BETA_F) / np.sqrt(SE_M**2 + SE_F**2 - 2 * r * SE_M * SE_F)

    # pvalue
    p = 2 * norm.sf(np.abs(t))

    return p


def getParser():
    parser = argparse.ArgumentParser(
        prog=str(__file__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
        @Author: xutingfeng@big.ac.cn 
        @Description: Run post-GWAS analysis

        Example:
        python3 P2.3-postGWAS.py -f GWASResultV2/Disease/Format/Gout -o GWASResultV2/Disease/PostGWAS 

        this will save to GWASResultV2/Disease/PostGWAS/Gout with all results 
        Contained steps:
        1. load sumstats from GWASResultV2/Disease
        2. assign counts to annotation
        3. genomic control
        4. manhattan plot
        5. trumpet plot
        6. miami plot
        7. heterogeneity test
        8. compare heterogeneity
        9. estimate heritability
        10. estimate genetic correlation

        """
        ),
    )
    # main params
    parser.add_argument(
        "-f",
        "--rootDir",
        type=str,
        help="rootDir",
        default="GWASResultV2/Disease/Format/Gout",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outputRottFolder",
        type=str,
        help="will save to thisFolder/rootDir.stem",
        default="GWASResultV2/Disease/PostGWAS",
        required=True,
    )
    parser.add_argument(
        "--ref_ld_chr",
        type=str,
        help="ref_ld_chr, e.g.:/pmaster/xutingfeng/dataset/ukb/dataset/LD_reference/LDSC/eur_w_ld_chr/ ; must end with / ",
        default="ref_ld_chr/",
        required=False,
    )
    parser.add_argument(
        "--female_counts_json_dir",
        type=str,
        help="female_counts_json_dir",
        default="white_female_counts.json",
        required=False,
    )
    parser.add_argument(
        "--male_counts_json_dir",
        type=str,
        help="male_counts_json_dir",
        default="white_male_counts.json",
        required=False,
    )
    parser.add_argument("--verbose", action="store_true", help="verbose", default=False)

    return parser


if __name__ == "__main__":
    log = Log()
    parser = getParser()
    args = parser.parse_args()
    verbose = args.verbose

    default_dpi = 400

    # config font
    import matplotlib.font_manager as fm
    import matplotlib.pyplot as plt

    # 添加自定义字体路径
    custom_font_path = "/home/xutingfeng/.fonts/Arial.ttf"  # 替换为你字体文件所在的目录
    fm.fontManager.addfont(custom_font_path)

    # check

    # 打印 Matplotlib 可识别的字体名称
    try:
        font_list = [f.name for f in fm.fontManager.ttflist]
        if "Arial" in font_list:
            # print("Arial font is available!")
            log.write("- Arial font is available!...", verbose=verbose)

        else:
            # print("Arial font is not found.")
            log.write("- Arial font is not found.", verbose=verbose)
    except:
        log.write("- Error in checking font", verbose=verbose)

    # step1: load
    rootDir = Path(args.rootDir)
    ref_ld_chr = args.ref_ld_chr
    trait = rootDir.stem
    saveDir = Path(args.outputRottFolder) / trait
    saveDir.mkdir(parents=True, exist_ok=True)

    getTraitDir = lambda sex: rootDir / f"{trait}_{sex}.feather"

    ## step2 assign conuts to annotation (required)
    # 从 JSON 文件读取统计结果
    female_counts_json_dir = args.female_counts_json_dir
    male_counts_json_dir = args.male_counts_json_dir

    if not Path(female_counts_json_dir).exists():
        raise FileNotFoundError(f"-Error with loading {female_counts_json_dir}")

    if not Path(male_counts_json_dir).exists():
        raise FileNotFoundError(f"-Error with loading {male_counts_json_dir}")

    # step2: gl.Sumstats
    try:
        Female_mysumstats = gl.Sumstats(
            pd.read_feather(getTraitDir("female")),
            snpid="SNPID",
            chrom="CHROM",
            pos="GENPOS",
            ea="ALLELE1",
            nea="ALLELE0",
            eaf="A1FREQ",
            n="N",
            # p="P",
            mlog10p="LOG10P",
            beta="BETA",
            se="SE",
            build="38",
            rsid="ID",
        )

        Male_mysumstats = gl.Sumstats(
            pd.read_feather(getTraitDir("male")),
            snpid="SNPID",
            chrom="CHROM",
            pos="GENPOS",
            ea="ALLELE1",
            nea="ALLELE0",
            eaf="A1FREQ",
            n="N",
            # p="P",
            mlog10p="LOG10P",
            beta="BETA",
            se="SE",
            rsid="ID",
            build="38",
        )
        ## load counts
        female_counts_json = json.load(open(female_counts_json_dir, "r"))
        male_counts_json = json.load(open(male_counts_json_dir, "r"))

        if trait in female_counts_json.keys() and trait in male_counts_json.keys():
            trait_type = "bt"
            # male

            n_case_male = male_counts_json[trait]["Cases"]
            n_control_male = male_counts_json[trait]["Controls"]
            # female
            n_case_female = female_counts_json[trait]["Cases"]
            n_control_female = female_counts_json[trait]["Controls"]
            count_dict = {
                "Female": {
                    "Case": n_case_female,
                    "Control": n_control_female,
                },
                "Male": {
                    "Case": n_case_male,
                    "Control": n_control_male,
                },
            }
        else:
            trait_type = "qt"
            count_dict = {
                "Female": {"N": Female_mysumstats.data["N"].iloc[0]},
                "Male": {"N": Male_mysumstats.data["N"].iloc[0]},
            }

        log.write(
            f"Starting {trait} analysis with trait type {trait_type}", verbose=verbose
        )
        log.write(f"count_dict is {count_dict}", verbose=verbose)

    except Exception as e:
        Log.write(f"Error with loading {trait} sumstats", verbose=verbose)
        Log.write(e, verbose=verbose)
        raise e

    # optional, but should be done there or before
    try:
        Female_mysumstats.remove_dup(mode="md")
        Male_mysumstats.remove_dup(mode="md")

        Female_mysumstats.data["P"] = 10 ** (-Female_mysumstats.data["MLOG10P"])
        Male_mysumstats.data["P"] = 10 ** (-Male_mysumstats.data["MLOG10P"])

        Female_mysumstats.data["Z"] = (
            Female_mysumstats.data["BETA"] / Female_mysumstats.data["SE"]
        )
        Male_mysumstats.data["Z"] = (
            Male_mysumstats.data["BETA"] / Male_mysumstats.data["SE"]
        )
        Female_mysumstats.data["MAF"] = Female_mysumstats.data["EAF"].apply(
            lambda x: min(x, 1 - x) if pd.notnull(x) else np.nan
        )
        Male_mysumstats.data["MAF"] = Male_mysumstats.data["EAF"].apply(
            lambda x: min(x, 1 - x) if pd.notnull(x) else np.nan
        )
    except Exception as e:
        log.write(f"Error with processing {trait} sumstats", verbose=verbose)
        raise e

    # step3: genomic control
    try:
        Female_lambda = calculate_lambda_from_p(Female_mysumstats.data["P"])
        Male_lambda = calculate_lambda_from_p(Male_mysumstats.data["P"])
        print(
            f"Female genomic inflation factor is {Female_lambda:.2f}\nMale genomic inflation factor is {Male_lambda:.2f}"
        )
        genomic_control = pd.DataFrame(
            {
                "Trait": [trait],
                "Female_genomic_inflation_factor": [Female_lambda],
                "Male_genomic_inflation_factor": [Male_lambda],
            }
        )

        genomic_control.to_csv(saveDir / "genomic_control.csv", index=False)
    except Exception as e:
        log.write(
            f"- Step3 Error with cal genomic control for {trait}", verbose=verbose
        )
        log.write(e, verbose=verbose)
        raise e

    # step4: manhattan plot of different sex
    try:
        title_female = f"{trait} Female"
        title_male = f"{trait} Male"

        if trait_type == "bt":
            title_female += f" (Case:{count_dict['Female']['Case']}, Control:{count_dict['Male']['Control']})"

            title_male += f" (Case:{count_dict['Female']['Case']}, Control:{count_dict['Male']['Control']})"

        else:
            title_male = f" Total N={count_dict['Male']['N']}"
            title_female = f" Total N={count_dict['Female']['N']}"

        Female_mahanttanPlot = Female_mysumstats.plot_mqq(
            skip=2, anno="GENENAME", title=title_female
        )
        Male_mahanttanPlot = Male_mysumstats.plot_mqq(
            skip=2, title=title_male, anno="GENENAME"
        )

        save_fig(
            fig=Female_mahanttanPlot[0],
            path=saveDir / "Female_mahanttanPlot",
            dpi=default_dpi,
        )
        save_fig(
            fig=Male_mahanttanPlot[0],
            path=saveDir / "Male_mahanttanPlot",
            dpi=default_dpi,
        )
    except Exception as e:
        log.write(f"- Step4 Error with cal manhattan plot for {trait}", verbose=verbose)
        log.write(e, verbose=verbose)
        raise e

    # step5: plot trumpetPlot for lead SNP
    try:
        sig_level = 5e-8
        Female_mysumstats_lead_data = Female_mysumstats.get_lead(
            sig_level=sig_level, gls=True, anno=True
        )
        Male_mysumstats_lead_data = Male_mysumstats.get_lead(
            sig_level=sig_level, gls=True, anno=True
        )

        Female_sig_loci_numbers = Female_mysumstats_lead_data.data.shape[0]
        Male_sig_loci_numbers = Male_mysumstats_lead_data.data.shape[0]
        if Female_sig_loci_numbers == 0:
            log.write(f"-Female {trait} no lead SNP", verbose=verbose)

        else:
            Female_mysumstats_lead_data.data.reset_index(drop=True).to_feather(
                saveDir / "Female_lead.feather"
            )
            log.write(
                f"-Female {trait} found {Female_sig_loci_numbers} loci with windowsizekb=500 and sig_level={sig_level:.2e}",
                verbose=verbose,
            )

        if Male_sig_loci_numbers == 0:
            log.write(f"-Male {trait} no lead SNP", verbose=verbose)
        else:
            Male_mysumstats_lead_data.data.reset_index(drop=True).to_feather(
                saveDir / "Male_lead.feather"
            )

            log.write(
                f"-Male {trait} found {Male_sig_loci_numbers} loci with windowsizekb=500 and sig_level={sig_level:.2e}",
                verbose=verbose,
            )

        if trait_type == "qt":
            female_trumpetPlot_args = {
                "sig_level": sig_level,
                "n": "N",
                "ts": [0.8],
                "mode": "q",
            }
            male_trumpetPlot_args = {
                "sig_level": sig_level,
                "n": "N",
                "ts": [0.8],
                "mode": "q",
            }
        elif trait_type == "bt":
            female_trumpetPlot_args = {
                "sig_level": sig_level,
                "ncase": count_dict["Female"]["Case"],
                "ncontrol": count_dict["Female"]["Control"],
                "prevalence": count_dict["Female"]["Case"] / count_dict["Female"]["N"],
                "ts": [0.8],
                "mode": "b",
            }
            male_trumpetPlot_args = {
                "sig_level": sig_level,
                "ncase": count_dict["Male"]["Case"],
                "ncontrol": count_dict["Male"]["Control"],
                "prevalence": count_dict["Male"]["Case"] / count_dict["Male"]["N"],
                "ts": [0.8],
                "mode": "b",
            }

        Female_mysumstats_lead_data_trumpetPlot = Female_mysumstats_lead_data.plot_trumpet(
            # sig_level=5e-6,
            # p_level=5e-6,
            anno="GENENAME",
            anno_style="right",
            cmap="tab20",
            or_to_rr=True,
            build="38",
            anno_x=0.01,
            anno_y=0,
            n_matrix=2000,
            fontsize=12,
            xscale="log",
            repel_force=0.15,
            sort="eaf",
            # hue="GENE",
            ylim=(-5, 4),
            **female_trumpetPlot_args,
            # title="Female",
        )
        Male_mysumstats_lead_data_trumpetPlot = Male_mysumstats_lead_data.plot_trumpet(
            # ts=[0.2, 0.4, 0.6, 0.8],
            anno="GENENAME",
            anno_style="right",
            cmap="tab20",
            or_to_rr=True,
            build="38",
            anno_x=0.01,
            anno_y=0,
            n_matrix=2000,
            fontsize=12,
            xscale="log",
            repel_force=0.15,
            sort="eaf",
            # hue="GENE",
            ylim=(-5, 4),
            **male_trumpetPlot_args,
        )

        save_fig(
            Male_mysumstats_lead_data_trumpetPlot,
            saveDir / "Male_trumpetPlot",
            dpi=default_dpi,
        )
        save_fig(
            Female_mysumstats_lead_data_trumpetPlot,
            saveDir / "Female_trumpetPlot",
            dpi=default_dpi,
        )
    except Exception as e:
        log.write(f"- Step5 Error with cal trumpet plot for {trait}", verbose=verbose)
        log.write(e, verbose=verbose)
        raise e

    # step6: miami plot
    try:
        Female_Male_miamiPlot = gl.plot_miami2(
            path1=Male_mysumstats,
            path2=Female_mysumstats,
            skip=2,
            #  id1="SNPID",
            #  id2="SNPID",
            anno1="GENENAME",
            anno2="GENENAME",
            build="38",
            titles=["Male" + title_male, "Female" + title_female],
        )

        save_fig(Female_Male_miamiPlot[0], saveDir / "miamiPlot", dpi=default_dpi)
    except Exception as e:
        log.write(f"-Step6 Error with cal miami plot for {trait}", verbose=verbose)
        log.write(e, verbose=verbose)
        raise e

    # step7: cal heterogeneity
    try:
        ## extract inner set
        innerSet = Female_mysumstats.data.merge(
            Male_mysumstats.data,
            on=["SNPID", "CHR", "POS", "EA", "NEA"],
            suffixes=["_female", "_male"],
        )

        ## cal spearman correlation
        from scipy.stats import spearmanr

        r, p = spearmanr(innerSet["P_male"], innerSet["P_female"])
        log.write(
            f"-Spearman rank correlation cross sex is {r:.4f}, p-value is {p:.4e}",
            verbose=verbose,
        )

        innerSet["het_p"] = het_test(
            BETA_F=innerSet["BETA_female"],
            SE_F=innerSet["SE_female"],
            BETA_M=innerSet["BETA_male"],
            SE_M=innerSet["SE_male"],
        )

        SexDiff_sumstats = gl.Sumstats(
            innerSet,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="het_p",
            build="38",
        )

        SexDiff_mahanttanPlot = SexDiff_sumstats.plot_mqq(
            skip=2, title=f"{trait} Sex Heterogeneity", anno="GENENAME"
        )

        save_fig(
            SexDiff_mahanttanPlot[0], saveDir / "SexDiff_mahanttanPlot", dpi=default_dpi
        )

        innerSet.reset_index(drop=True).to_feather(saveDir / "sex_het.feather")

        SigLoci = SexDiff_sumstats.get_lead(sig_level=1e-8, anno=True, gls=True)

        if SigLoci.shape[0] == 0:
            log.write(f"-No significant loci found of {trait}", verbose=verbose)
            SigLoci.data.to_feather("sex_het_SigLoci.feather")
        else:
            SigLoci.data.reset_index(drop=True).to_feather("sex_het_SigLoci.feather")
    except Exception as e:
        log.write(f"-Step7 Error with cal heterogeneity for {trait}", verbose=verbose)
        log.write(e, verbose=verbose)
        raise e

    # step8: compare het plot
    # try:
    #     female_male_comparePlot = gl.compare_effect(
    #         path1=Female_mysumstats,
    #         cols_name_list_1=["SNPID", "P", "NEA", "EA", "CHR", "POS"],
    #         effect_cols_list_1=["BETA", "SE"],
    #         path2=Male_mysumstats,
    #         cols_name_list_2=["SNPID", "P", "NEA", "EA", "CHR", "POS"],
    #         effect_cols_list_2=["BETA", "SE"],
    #         label=["Female", "Male", "Both", "None"],
    #         xylabel_prefix="Per-allele effect size for ",
    #         anno="GENENAME",
    #         anno_het=True,
    #         anno_diff=0.015,
    #         is_q=True,
    #         sig_level=5e-8,
    #         legend_title=r"$ P < 5 x 10^{-8}$ in:",
    #         verbose=True,
    #         build="38",
    #         # mode="OR",
    #         #   save = "myplot.png",
    #         #   save_args= {"dpi":300,"facecolor":"white"}
    #     )
    #     save_fig(
    #         female_male_comparePlot[0],
    #         saveDir / "female_male_comparePlot",
    #         dpi=default_dpi,
    #     )
    # except Exception as e:
    #     log.write(
    #         f"-Step8 Error with cal compare het plot for {trait}", verbose=verbose
    #     )
    #     log.write(e, verbose=verbose)
    #     raise e

    # step9 estimate heritability
    try:
        ref_ld_chr = ref_ld_chr + "/"  # to ensure end with /
        w_ld_chr = ref_ld_chr

        # keep only hapmap3 SNPs

        Female_mysumstats.filter_hapmap3(inplace=True)
        Male_mysumstats.filter_hapmap3(inplace=True)

        Female_mysumstats.data.reset_index(drop=True).to_feather(
            saveDir / "Female_hm3.feather"
        )
        Male_mysumstats.data.reset_index(drop=True).to_feather(
            saveDir / "Male_hm3.feather"
        )

        Female_mysumstats.estimate_h2_by_ldsc(ref_ld_chr=ref_ld_chr, w_ld_chr=w_ld_chr)
        Male_mysumstats.estimate_h2_by_ldsc(ref_ld_chr=ref_ld_chr, w_ld_chr=w_ld_chr)

        h2_df = pd.concat(
            [
                Male_mysumstats.ldsc_h2.assign(Trait=trait, Group="Male").drop(
                    columns=["Catagories"]
                ),
                Female_mysumstats.ldsc_h2.assign(Trait=trait, Group="Female").drop(
                    columns=["Catagories"]
                ),
            ]
        )
        # format the dataframe
        h2_df = h2_df.pivot(index="Trait", columns="Group")
        # flatten the column
        h2_df.columns = h2_df.columns.to_flat_index()
        h2_df.columns = ["_".join(col).strip() for col in h2_df.columns.values]
        # reset index
        h2_df = h2_df.reset_index()
        # set dtype
        # h2_df.iloc[:, 1:] = h2_df.iloc[:, 1:].astype(float)
        for col in h2_df.columns[1:]:
            h2_df[col] = h2_df[col].astype(float)

        # set the heterogeneity test; note this is a two-tailed test with largely sample size, so think it is normal distribution
        h2_df["het_p"] = het_test(
            BETA_F=h2_df["h2_obs_Female"],
            SE_F=h2_df["h2_se_Female"],
            BETA_M=h2_df["h2_obs_Male"],
            SE_M=h2_df["h2_se_Male"],
            r=0,
        )

        h2_df.to_csv(saveDir / "h2_df.csv", index=False)
    except Exception as e:
        log.write(f"-Step9 Error with cal heritability for {trait}", verbose=verbose)
        log.write(e, verbose=verbose)
        raise e

    # step10: cal genetic correlation
    try:

        Female_mysumstats.estimate_rg_by_ldsc(
            other_traits=[Male_mysumstats],
            rg="Female,Male",
            ref_ld_chr=ref_ld_chr,
            w_ld_chr=w_ld_chr,
        )
        ldsc_rg_df = Female_mysumstats.ldsc_rg.copy()
        ldsc_rg_df.insert(0, "Trait", trait)

        rg = Female_mysumstats.ldsc_rg["rg"][0]
        se_rg = Female_mysumstats.ldsc_rg["se"][0]

        t_sex_diff = (rg - 1) / se_rg

        ldsc_rg_df["p_rg_diff"] = cal_p_from_z(t_sex_diff)
        ldsc_rg_df.to_csv(saveDir / "ldsc_rg_df.csv", index=False)

        log.write(f"-Ending {trait} analysis", verbose=verbose)
    except Exception as e:
        log.write(
            f"-Step10 Error with cal genetic correlation for {trait}", verbose=verbose
        )
        log.write(e, verbose=verbose)
        raise e
