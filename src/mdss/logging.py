import logging
import os
import sys
import os.path
from datetime import datetime

here = os.path.abspath(os.path.dirname(__file__))
log_dir = os.path.join(here, "logs")

logging_configured = False

if not logging_configured:
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    filename = datetime.now().strftime("log_%Y_%m_%d_%H_%M_%S.txt")
    filepath = os.path.join(log_dir, filename)

    logging.basicConfig(
        filename=filepath,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        filemode="w",
    )

    logging_configured = True
