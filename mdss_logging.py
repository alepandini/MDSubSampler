import logging

logging.basicConfig(
    filename="log.txt",
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    filemode="w",
)
logging.debug("Logging test...")
logging.info("The program is working as expected")
logging.warning("The program may not function properly")
logging.error("The program encountered an error")
logging.critical("The program crashed")
