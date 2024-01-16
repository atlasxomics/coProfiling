import os
import glob
import subprocess
from enum import Enum
from pathlib import Path
from dataclasses import dataclass
from typing import List, Optional, Union, Tuple
from latch.registry.table import Table
from latch.resources.launch_plan import LaunchPlan
from latch import large_task, large_task, small_task, workflow

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
    gex_dir: LatchDir = LatchDir("latch:///star_outputs/demo/STAR_outputsGene/raw")
    fragments_file: LatchFile = LatchFile("latch:///chromap_outputs/demo/chromap_output/fragments.tsv.gz")
    spatial_dir: LatchDir = LatchDir("latch:///spatials/demo/spatial/")
    
    

class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'

@large_task(retries=0)
def coProf_task(runs: List[Run],
                project_name: str,
                genome: Genome,
                spot_size: int,
                cluster_resolution: float
                ) -> LatchDir:
    output_dir = Path("reports/").resolve()
    os.mkdir(output_dir)
    _report_cmd = [
            "Rscript",
            "/root/wf/coProf.R",
            project_name,
            genome.value,
            f'{spot_size}',
            f'{cluster_resolution}'
        
    ]
    runs = [
        (
            f'{run.run_id},'
            f'{run.gex_dir.local_path},'
            f'{run.fragments_file.local_path},'
            f'{run.spatial_dir.local_path},'
        )
        for run in runs
    ]

    _report_cmd.extend(runs)
    subprocess.run(_report_cmd)
    return LatchDir(str(output_dir), f"latch:///coProfiling_outputs/{project_name}/reports")

metadata = LatchMetadata(
    display_name='coProfiling_wf',
    author=LatchAuthor(
        name='AtlasXomics, Inc.',
        email='nouroddinc@atlasxomics.com',
        github='github.com/atlasxomics/coProfiling',
    ),
    repository='https://github.com/atlasxomics/coProfiling',
    license='MIT',
    parameters={
        'runs': LatchParameter(
            display_name='runs',
            description='List of runs to be analyzed; each run must contain a \
                        run_id and fragments.tsv.gz file; gene expression and \
                        spatial folder for SpatialDimPlot.',
            batch_table_column=True,
            samplesheet=True
        ),
            'project_name': LatchParameter(
            display_name='project name',
            description='Name of output directory in coProfiling_outputs/',
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="project name cannot start with a '/'"
                )
            ]
        ),
        'genome': LatchParameter(
            display_name='genome',
            description='Reference genome to be used for geneAnnotation',
            batch_table_column=True,
        ),
        'spot_size': LatchParameter(
            display_name='spot size',
            description='The size of the spots in the spatial maps',
            batch_table_column=True,
            hidden=True
        ),
        'cluster_resolution': LatchParameter(
            display_name='cluster resolution',
            description='resolution parameter for FindClusters function',
            batch_table_column=True,
            hidden=True
        )},
)
@workflow(metadata)
def coProfiling_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    spot_size: int = 1,
    cluster_resolution: float = 0.5
) -> LatchDir:
        reports = coProf_task(
        runs=runs,
        project_name=project_name,
        genome=genome,
        spot_size=spot_size,
        cluster_resolution=cluster_resolution,

    )
        return reports

LaunchPlan(
    coProfiling_workflow,
    'defaults',
    {
        'runs': [
            Run(
                'default',
                LatchDir('latch:///star_outputs/demo/STAR_outputsGene/raw'),
                LatchFile('latch:///chromap_outputs/demo/chromap_output/fragments.tsv.gz'),
                LatchDir('latch:///spatials/demo/spatial/')
                )
        ],
        'project_name': 'coProf',
        'genome': Genome.hg38,
        'spot_size': 1,
        'cluster_resolution': 0.5
    },
)


