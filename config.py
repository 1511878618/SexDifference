from pathlib import Path 
version = "v1"

# Path to the project root directory

projectDir = Path(version)
projectDir.mkdir(exist_ok=True, parents=True)
# for figure
paperDir = projectDir / "paper"
paperDir.mkdir(exist_ok=True)
paperFigDir = paperDir / "fig"
paperFigDir.mkdir(exist_ok=True)
paperTableDir = paperDir / "table"
paperTableDir.mkdir(exist_ok=True)

# raw data 
dataDir = projectDir / "data"
dataDir.mkdir(exist_ok=True)
outputDir = projectDir / "output"
outputDir.mkdir(exist_ok=True)

