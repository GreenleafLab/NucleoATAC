#Define version based on setup script
import pkg_resources
__version__ = pkg_resources.require("NucleoATAC")[0].version
