import os
import subprocess

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List

from latch.resources.launch_plan import LaunchPlan
from latch import large_task, workflow

from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchFile,
    LatchMetadata,
    LatchParameter,
    LatchRule
)


@dataclass
class Run:
    run_id: str
    gex_dir: LatchDir = LatchDir(
        "latch:///star_outputs/demo/STAR_outputsGene/raw"
    )
    fragments_file: LatchFile = LatchFile(
        "latch:///chromap_outputs/demo/chromap_output/fragments.tsv.gz"
    )
    spatial_dir: LatchDir = LatchDir(
        "latch:///spatials/demo/spatial/"
    )


class Genome(Enum):
    mm10 = "mm10"
    hg38 = "hg38"


@large_task(retries=0)
def coProf_task(
    runs: List[Run],
    project_name: str,
    genome: Genome,
    spot_size: int,
    cluster_resolution: float,
    tile_size: int,
    min_TSS: float,
    min_frags: int,
    lsi_iterations: int,
    lsi_resolution: float,
    lsi_varfeatures: int
) -> LatchDir:

    output_dir = Path("reports/").resolve()
    os.makedirs(output_dir, exist_ok=True)

    # Turn on virtual dev for png in Rmd
    _xvfb_cmd = ["Xvfb", ":99", "-screen", "0", "1024x768x16"]
    subprocess.Popen(_xvfb_cmd)
    os.environ["DISPLAY"] = ":99"

    _report_cmd = [
        "Rscript",
        "/root/wf/coProf.R",
        project_name,
        genome.value,
        f"{spot_size}",
        f"{cluster_resolution}",
        f"{tile_size}",
        f"{min_TSS}",
        f"{min_frags}",
        f"{lsi_iterations}",
        f"{lsi_resolution}",
        f"{lsi_varfeatures}"
    ]
    runs = [
        (
            f"{run.run_id},"
            f"{run.gex_dir.local_path},"
            f"{run.fragments_file.local_path},"
            f"{run.spatial_dir.local_path},"
        )
        for run in runs
    ]

    _report_cmd.extend(runs)
    subprocess.run(_report_cmd)

    return LatchDir(
        str(output_dir),
        f"latch:///coProfiling_outputs/{project_name}/reports"
    )


metadata = LatchMetadata(
    display_name="coProfiling_wf",
    author=LatchAuthor(
        name="AtlasXomics, Inc.",
        email="nouroddinc@atlasxomics.com",
        github="github.com/atlasxomics/coProfiling",
    ),
    repository="https://github.com/atlasxomics/coProfiling",
    license="MIT",
    parameters={
        "runs": LatchParameter(
            display_name="runs",
            description="List of runs to be analyzed; each run must contain a \
                        run_id and fragments.tsv.gz file; gene expression \
                        and spatial folder for SpatialDimPlot.",
            batch_table_column=True,
            samplesheet=True
        ),
        "project_name": LatchParameter(
            display_name="project name",
            description="Name of output directory in coProfiling_outputs/",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="project name cannot start with a '/'"
                )
            ]
        ),
        "genome": LatchParameter(
            display_name="genome",
            description="Reference genome to be used for geneAnnotation",
            batch_table_column=True,
        ),
        "spot_size": LatchParameter(
            display_name="spot size",
            description="The size of the spots in the spatial maps",
            batch_table_column=True,
            hidden=True
        ),
        "cluster_resolution": LatchParameter(
            display_name="cluster resolution",
            description="resolution parameter for FindClusters function",
            batch_table_column=True,
            hidden=True
        ),
        "tile_size": LatchParameter(
            display_name="tile size",
            description="The size of the tiles used for binning counts in the \
                        TileMatrix.",
            batch_table_column=True,
            hidden=True
        ),
        "min_TSS": LatchParameter(
            display_name="minimum TSS",
            description="The minimum numeric transcription start site (TSS) \
                        enrichment score required for a cell to pass \
                        filtering.",
            batch_table_column=True,
            hidden=True
        ),
        "min_frags": LatchParameter(
            display_name="minimum fragments",
            description="The minimum number of mapped ATAC-seq fragments \
                        required per cell to pass filtering.",
            batch_table_column=True,
            hidden=True
        ),
        "lsi_iterations": LatchParameter(
            display_name="LSI iterations",
            description="iterations parameter from addIterativeLSI function.",
            batch_table_column=True,
            hidden=True
        ),
        "lsi_resolution": LatchParameter(
            display_name="LSI resolution",
            description="resolution parameter from \
                        addIterativeLSI/clusterParams function.",
            batch_table_column=True
        ),
        "lsi_varfeatures": LatchParameter(
            display_name="LSI varFeatures",
            description="varFeatures parameter from addIterativeLSI function; \
                        each will correspond to a umap.pdf, the last in the \
                        will be used to make the RDS object.",
            batch_table_column=True
        )},
)


@workflow(metadata)
def coProfiling_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    spot_size: int = 1,
    cluster_resolution: float = 0.5,
    tile_size: int = 5000,
    min_TSS: float = 2.0,
    min_frags: int = 0,
    lsi_iterations: int = 2,
    lsi_resolution: float = 0.5,
    lsi_varfeatures: int = 25000
) -> LatchDir:

    reports = coProf_task(
        runs=runs,
        project_name=project_name,
        genome=genome,
        spot_size=spot_size,
        cluster_resolution=cluster_resolution,
        tile_size=tile_size,
        min_TSS=min_TSS,
        min_frags=min_frags,
        lsi_iterations=lsi_iterations,
        lsi_resolution=lsi_resolution,
        lsi_varfeatures=lsi_varfeatures
    )
    return reports


LaunchPlan(
    coProfiling_workflow,
    "defaults",
    {
        "runs": [
            Run(
                "default",
                LatchDir("latch:///star_outputs/demo/STAR_outputsGene/raw"),
                LatchFile("latch://13502.account/sample_fqs/test/fragments.tsv.gz"),
                LatchDir("latch://atx-illumina.mount/Images_spatial/D1372/spatial")
            )
        ],
        "project_name": "demo",
        "genome": Genome.mm10,
        "spot_size": 1,
        "cluster_resolution": 0.5,
        "tile_size": 5000,
        "min_TSS": 2.0,
        "min_frags": 0,
        "lsi_iterations": 2,
        "lsi_resolution": 0.5,
        "lsi_varfeatures": 25000
    },
)

if __name__ == "__main__":

    coProf_task(
        runs=[
            Run(
                run_id="D1372",
                gex_dir="latch://13502.account/star_outputs/merged_D1372_D1373/STAR_outputs/D01372_NG02965_STAR_outputsGene/raw",
                fragments_file="latch://13502.account/noori_sample_fqs/test/fragments.tsv.gz",
                spatial_dir="latch://atx-illumina.mount/Images_spatial/D1372/spatial"
            )
        ],
        project_name="D1372_dev",
        genome=Genome("mm10"),
        spot_size=1,
        cluster_resolution=0.5,
        tile_size=5000,
        min_TSS=2.0,
        min_frags=0,
        lsi_iterations=2,
        lsi_resolution=0.5,
        lsi_varfeatures=25000
    )
