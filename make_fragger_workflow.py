from pathlib import Path

def update_fasta_in_workflow(
    workflow_file: str,
    new_fasta: str,
    output_file: str | None = None
) -> str:
    """
    Replace database.db-path in a FragPipe workflow file.

    Parameters
    ----------
    workflow_file : str
        Existing workflow file.
    new_fasta : str
        Path to the new FASTA file.
    output_file : str | None
        Output workflow filename. If None, creates
        '<workflow_file>.updated'.

    Returns
    -------
    str
        Path to the written workflow file.
    """

    workflow_file = Path(workflow_file)

    if output_file is None:
        output_file = workflow_file.with_suffix(
            workflow_file.suffix + ".updated"
        )

    found = False
    output_lines = []

    with open(workflow_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("database.db-path="):
                output_lines.append(
                    f"database.db-path={new_fasta}\n"
                )
                found = True
            else:
                output_lines.append(line)

    if not found:
        raise ValueError(
            "database.db-path entry not found in workflow file"
        )

    with open(output_file, "w", encoding="utf-8") as f:
        f.writelines(output_lines)

    return str(output_file)

if __name__ == '__main__':
    import sys
    new_workflow = update_fasta_in_workflow(
    workflow_file=sys.argv[1],
    new_fasta=sys.argv[2],
    output_file=sys.argv[3])
