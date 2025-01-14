#!/usr/bin/env python


import pandas as pd
import subprocess
import time
import io
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import textwrap
import shutil


def check_command_exist(commands):
    if isinstance(commands, str):
        commands = [commands]

    errors = ""
    for cmd in commands:
        if not shutil.which(cmd):
            errors += f"{cmd} not found in $PATH\n"

    if len(errors) > 0:
        raise FileNotFoundError(errors)
    return True


def getParser():
    parser = argparse.ArgumentParser(
        prog=str(__file__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
        @Author: xutingfeng@big.ac.cn 
        @Description: Merge the regenie results of the same trait from different chromosomes


        Assume this for cmd code, required for resetID2.py and versionConvert.py for gwas 
        CHR POS ID A0(REF) A1(ALT)

        

        Example:
            # For GWAS Results
            python %(prog)s --tgtDir ./GWASResult/RF 
            # For Burden Results
            python %(prog)s --tgtDir ./BurdenResult/RF --is_burden

        """
        ),
    )
    # main params
    parser.add_argument(
        "--tgtDir",
        type=str,
        help="The target trait directory, this should contain two subfolders, 1) female and 2) male, each contains the regenie results of the trait for different chromosomes, tgtDir.name should be the trait name or passed by --tgt",
    )
    parser.add_argument(
        "--is_burden",
        action="store_true",
        help="Whether the burden test",
    )
    parser.add_argument(
        "--saveDir",
        type=str,
        default=None,
        help="The save directory, if not provided, the results will be saved in the tgtDir/merged",
    )

    parser.add_argument(
        "-c", "--chain", default=[], nargs="+", help="Chain file, like -c hg19 hg38"
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=4, help="The number of threads to use"
    )

    return parser


# 定义并行执行的任务函数，接收字典参数
def process_trait(trait_params):
    trait = trait_params["trait"]
    trait_list_df = trait_params["trait_list_df"]
    saveDir = trait_params["saveDir"]
    versionConvert = trait_params.get("versionConvert", False)
    is_burden = trait_params.get("is_burden", False)
    from_version = trait_params.get("from_version", None)
    to_version = trait_params.get("to_version", None)

    for sex, chr_list_df in trait_list_df.groupby("Sex"):
        trait_save_folder = saveDir / trait
        trait_save_folder.mkdir(parents=True, exist_ok=True)

        sex_trait_results_saveDir = trait_save_folder / f"{trait}_{sex}.feather"

        time_start = time.time()
        print(f"Processing {sex} {trait}")

        if sex_trait_results_saveDir.exists():
            continue
        sex_res_list = []
        for regenie_result_file in tqdm(
            chr_list_df["file_dir"], desc=f"Loading {sex} {trait}"
        ):
            #
            if not is_burden:

                cmd = f"cat {regenie_result_file}"

                if versionConvert:
                    cmd += (
                        f"|versionConvert.py -c {from_version} {to_version} -i 1 2 -n --drop"
                    )

                cmd += "| resetID2.py --add-col SNPID -i 1 2 5 4"
                # print(cmd)
                result = subprocess.run(
                    cmd,
                    text=True,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )

                if result.returncode != 0:
                    print(f"Command failed with return code {result.returncode}")
                    print(f"Error output: {result.stderr}")
                    raise RuntimeError("Pipeline command failed.")

                trait_chr_res_df = pd.read_csv(io.StringIO(result.stdout), sep=r"\s+")
                print(f"Succeed in converting {regenie_result_file}")
            else:
                trait_chr_res_df = pd.read_csv(regenie_result_file, sep=r"\s+")

            sex_res_list.append(trait_chr_res_df)

        trait_res_df = (
            pd.concat(sex_res_list).drop(columns=["EXTRA"]).reset_index(drop=True)
        )

        # trait_res_df.to_parquet(sex_trait_results_saveDir)
        import pdb; pdb.set_trace()
        trait_res_df.to_feather(sex_trait_results_saveDir)
        time_end = time.time()
        # Duration in seconds
        duration = time_end - time_start
        print(
            f"Finished {sex} {trait} and saved to {sex_trait_results_saveDir} with time {duration:.2f} s"
        )


# 并行化执行
def parallel_process_traits(
    files_df,
    saveDir,
    is_burden=False,
    max_workers=4,
    versionConvert=False,
    from_version=None,
    to_version=None,
):
    if max_workers == 1:
        for trait, trait_list_df in files_df.groupby("Name"):
            process_trait(
                {
                    "trait": trait,
                    "trait_list_df": trait_list_df,
                    "saveDir": saveDir,
                    "is_burden": is_burden,
                    "versionConvert": versionConvert,
                    "from_version": from_version,
                    "to_version": to_version,
                }
            )
    else:
        # 为每个 trait 创建参数字典
        traits = [
            {
                "trait": trait,
                "trait_list_df": trait_list_df,
                "saveDir": saveDir,
                "is_burden": is_burden,
                "versionConvert": versionConvert,
                "from_version": from_version,
                "to_version": to_version,
            }
            for trait, trait_list_df in files_df.groupby("Name")
        ]

        # 使用 tqdm.rich 显示进度条
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # 提交任务并显示进度条
            list(
                tqdm(
                    executor.map(process_trait, traits),
                    total=len(traits),
                    desc="Processing Traits",
                )
            )


# params
parser = getParser()
args = parser.parse_args()
tgt_dir = Path(args.tgtDir)
is_burden = args.is_burden
saveDir = args.saveDir
threads = args.threads
# hg19ToHg38 = args.hg19ToHg38
chains = args.chain

# check
if len(chains) == 2:
    check_command_exist(["resetID2.py", "versionConvert.py"])
    versionConvert = True
    from_version, to_version = chains
elif len(chains) != 0 or len(chains) != 2:
    raise ValueError(
        "Chain file should be provided in pair, like -c hg19 hg38, but now is:",
        " ".join(chains),
    )


# save Dir
saveDir = tgt_dir / "merged" if saveDir is None else Path(saveDir)
saveDir.mkdir(parents=True, exist_ok=True)

if is_burden:
    combination_index = ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "TEST"]
else:
    combination_index = ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1"]


## 2)Load the results
files_list = []
for sex in ["female", "male"]:

    sex_result_dir = tgt_dir / sex

    regenie_results_files = [f for f in sex_result_dir.glob("*.regenie") if f.is_file()]

    for regenie_result_file in regenie_results_files:
        stem = regenie_result_file.stem
        chr_name = stem.split("_")[0]
        full_name = "_".join(stem.split("_")[1:])
        files_list.append(
            {
                "Name": full_name,
                "Chr": chr_name,
                "Sex": sex,
                "file_dir": str(regenie_result_file),
            }
        )

files_df = pd.DataFrame(files_list)

parallel_process_traits(
    files_df,
    saveDir,
    is_burden=is_burden,
    versionConvert=versionConvert,
    from_version=from_version,
    to_version=to_version,
    max_workers=threads,
)
