import hashlib
import shlex
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import click
from loguru import logger

"""
to make the BLAST dbs

# v 2.14.0

 makeblastdb -in dnaA.faa -dbtype prot -out dnaA_db
 makeblastdb -in terL.faa -dbtype prot -out terL_db
 makeblastdb -in repA.faa -dbtype prot -out repA_db
 makeblastdb -in all.faa -dbtype prot -out all_db

 makeblastdb -in dnaA_repA.faa -dbtype prot -out dnaA_repA_db
 makeblastdb -in dnaA_terL.faa -dbtype prot -out dnaA_terL_db
 makeblastdb -in repA_terL.faa -dbtype prot -out repA_terL_db


"""


class ExternalTool:
    def __init__(self, tool: str, input: str, output: str, params: str, logdir: Path):
        self.command: List[str] = self._build_command(tool, input, output, params)
        logdir.mkdir(parents=True, exist_ok=True)
        command_hash = hashlib.sha256(self.command_as_str.encode("utf-8")).hexdigest()
        tool_name = Path(tool).name
        logfile_prefix: Path = logdir / f"{tool_name}_{command_hash}"
        self.out_log = f"{logfile_prefix}.out"
        self.err_log = f"{logfile_prefix}.err"

    @property
    def command_as_str(self) -> str:
        return shlex.join(self.command)

    @staticmethod
    def _build_command(tool: str, input: str, output: str, params: str) -> List[str]:
        command = f"{tool} {input} {output} {params}"  # this is how mmseqs does it
        print(command)
        escaped_command = shlex.split(command)
        return escaped_command

    def run(self) -> None:
        with open(self.out_log, "w") as stdout_fh, open(self.err_log, "w") as stderr_fh:
            logger.info(f"Started running {self.command_as_str} ...")
            self._run_core(self.command, stdout_fh=stdout_fh, stderr_fh=stderr_fh)
            logger.info(f"Done running {self.command_as_str}")

    @staticmethod
    def _run_core(command: List[str], stdout_fh, stderr_fh) -> None:
        subprocess.check_call(command, stdout=stdout_fh, stderr=stderr_fh)

    @staticmethod
    def run_tools(
        tools_to_run: Tuple["ExternalTool", ...], ctx: Optional[click.Context] = None
    ) -> None:
        for tool in tools_to_run:
            try:
                tool.run()
            except subprocess.CalledProcessError as error:
                logger.error(
                    f"Error calling {tool.command_as_str} (return code {error.returncode})"
                )
                logger.error(f"Please check stdout log file: {tool.out_log}")
                logger.error(f"Please check stderr log file: {tool.err_log}")
                logger.error("Temporary files are preserved for debugging")
                logger.error("Exiting...")

                if ctx:
                    ctx.exit(1)
                else:
                    sys.exit(1)

    """
    Only one toolf
    """

    @staticmethod
    def run_tool(tool: "ExternalTool", ctx: Optional[click.Context] = None) -> None:
        try:
            tool.run()
        except subprocess.CalledProcessError as error:
            logger.error(
                f"Error calling {tool.command_as_str} (return code {error.returncode})"
            )
            logger.error(f"Please check stdout log file: {tool.out_log}")
            logger.error(f"Please check stderr log file: {tool.err_log}")
            logger.error("Temporary files are preserved for debugging")
            logger.error("Exiting...")

            if ctx:
                ctx.exit(1)
            else:
                sys.exit(1)
