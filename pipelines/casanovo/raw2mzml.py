import subprocess 
from pathlib import Path
import argparse


"""
Usage example: python raw2mzml.py -i BBM_428_P110_31_MIA_007_30.raw -o casanovo_input

"""


def main():
    parser = argparse.ArgumentParser(
        description="Convert Thermo RAW files to mzML using msconvert in Docker."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input RAW file"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory"
    )

    parser.add_argument(
        "--docker-image",
        default="chambm/pwiz-skyline-i-agree-to-the-vendor-licenses",
        help="Docker image containing msconvert"
    )

    args = parser.parse_args()

    run_msconvert_docker(
        rawfile_path=args.input,
        output_dir=args.output,
        docker_image=args.docker_image,
    )

def run_msconvert_docker(
    rawfile_path: str | Path,
    output_dir: str | Path,
    docker_image: str = "chambm/pwiz-skyline-i-agree-to-the-vendor-licenses",
) -> None:
    """
    Run ProteoWizard msconvert in Docker on Linux.
    """

    rawfile_path = Path(rawfile_path).resolve()
    output_dir = Path(output_dir).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)

    title_maker = (
        'titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState> '
        'File:"^<SourcePath^>", NativeID:"^<Id^>"'
    )

    # Mount the RAW file directory
    raw_dir = rawfile_path.parent

    cmd = [
        "docker",
        "run",
        "--rm",

        # RAW file directory
        "-v",
        f"{raw_dir}:/input",

        # Output directory
        "-v",
        f"{output_dir}:/output",

        docker_image,

        "wine",
        "msconvert",

        f"/input/{rawfile_path.name}",

        "--zlib",
        "--simAsSpectra",

        "--filter",
        "peakPicking vendor msLevel=1-",

        "--filter",
        title_maker,

        "--outdir",
        "/output",
    ]


    subprocess.run(cmd, check=True)





if __name__ == '__main__':
    main()
