"""
    @release_date  : $release_date
    @version       : $release_version
    @author        : Namir Oues
    
    This file is part of the MDSubSampler software 
    (https://github.com/alepandini/MDSubSampler).
    Copyright (c) 2023 Namir Oues and Alessandro Pandini.

    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
    General Public License for more details.

    You should have received a copy of the GNU General Public License 
    along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
import logging as log
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

    log.basicConfig(
        filename=filepath,
        level=log.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        filemode="w",
    )

    logging_configured = True
