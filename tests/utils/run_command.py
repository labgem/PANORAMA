import logging
import subprocess

# create a module-level logger
logger = logging.getLogger(__name__)


# def run_command(cmd):
#     logger.info(cmd)
#     with patch("sys.argv", cmd.split()):
#         main()


def run_command(cmd):
    logger.info(cmd)

    subprocess.run(cmd.split())
