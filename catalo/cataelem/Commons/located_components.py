# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

from cataelem.Tools.base_objects import LocatedComponents, objects_from_context
import cataelem.Commons.physical_quantities as PHY

# ----------------------------------------------------------------------------------------------
# Located components - ELEM - Field on element (constant)
# ----------------------------------------------------------------------------------------------
# Field for geometry (2D)
CGEOM2D = LocatedComponents(phys=PHY.GEOM_R, type="ELEM", components=("X", "Y"))

# Field for geometry (3D)
CGEOM3D = LocatedComponents(phys=PHY.GEOM_R, type="ELEM", components=("X", "Y", "Z"))

# Field for displacements (3D)
CDEPL3D = LocatedComponents(phys=PHY.DEPL_R, type="ELEM", components=("DX", "DY", "DZ"))

# Field for strains (Function-3D)
CEPS3DF = LocatedComponents(
    phys=PHY.EPSI_F, type="ELEM", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)

# Field for strains (Real-3D)
CEPS3DR = LocatedComponents(
    phys=PHY.EPSI_R, type="ELEM", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)

# Field for strains (Real-2D)
CEPS2DR = LocatedComponents(
    phys=PHY.EPSI_R, type="ELEM", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)

# Field for strains (Complex-2D)
CEPS2DC = LocatedComponents(
    phys=PHY.EPSI_C, type="ELEM", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)

# Field for strains (Function-2D)
CEPS2DF = LocatedComponents(
    phys=PHY.EPSI_F, type="ELEM", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)

# Field for distributed loads (function-2D)
CFOR2DF = LocatedComponents(phys=PHY.FORC_F, type="ELEM", components=("FX", "FY"))

# Field for distributed loads (function-3D)
CFOR3DF = LocatedComponents(phys=PHY.FORC_F, type="ELEM", components=("FX", "FY", "FZ"))

# Field for distributed loads (function-pipe)
CFORPIF = LocatedComponents(
    phys=PHY.FORC_F, type="ELEM", components=("FX", "FY", "FZ", "MX", "MY", "MZ", "REP")
)

# Field for distributed loads (Real-2D)
CFOR2DR = LocatedComponents(phys=PHY.FORC_R, type="ELEM", components=("FX", "FY"))

# Field for distributed loads (Real-3D)
CFOR3DR = LocatedComponents(phys=PHY.FORC_R, type="ELEM", components=("FX", "FY", "FZ"))

# Field for distributed loads (real-pipe)
CFORPIR = LocatedComponents(
    phys=PHY.FORC_R, type="ELEM", components=("FX", "FY", "FZ", "MX", "MY", "MZ", "REP")
)

# Field for pressure (Function-2D)
CPRE2DF = LocatedComponents(phys=PHY.PRES_F, type="ELEM", components=("PRES", "CISA"))

# Field for pressure (Function-3D)
CPRE3DF = LocatedComponents(phys=PHY.PRES_F, type="ELEM", components=("PRES",))

# Field for pressure (Real-3D)
CPRE3DR = LocatedComponents(phys=PHY.PRES_R, type="ELEM", components=("PRES",))

# Field for pressure (Real-3D for solid shell elements)
CPRESBR = LocatedComponents(phys=PHY.PRES_R, type="ELEM", components=("PINF", "PSUP"))

# Field for pressure (Function-3D)
CPRESBF = LocatedComponents(phys=PHY.PRES_F, type="ELEM", components=("PINF", "PSUP"))

# Field for 'EFFE_FOND'
CEFOND = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X1",))

# Output field for theta (fracture mechanic-XFEM-2D)
CKGTX2D = LocatedComponents(phys=PHY.RUPT_R, type="ELEM", components=("GTHETA", "FIC[2]", "K[2]"))

# Output field for theta (fracture mechanic-XFEM-3D)
CKGTX3D = LocatedComponents(
    phys=PHY.RUPT_R, type="ELEM", components=("GTHETA", "FIC[3]", "K[3]", "BETA")
)

# Output field for G-theta and K (fracture mechanic)
CKGTHET = LocatedComponents(phys=PHY.RUPT_R, type="ELEM", components=("GTHETA", "FIC[3]", "K[3]"))

# Output field for G-theta (fracture mechanic)
CGTHETA = LocatedComponents(phys=PHY.RUPT_R, type="ELEM", components=("GTHETA",))

# Field for abscisses of Lagrange polynom (fracture mechanic)
CABSLAG = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[3]",))

# Field for structural elements: radius
CRADIUS = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X1",))

# Field for structural elements: section (2D)
CSECT2D = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[6]",))

# Field for structural elements: section (3D)
CSECT3D = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[10]",))

# Field for L2-Norm
CNORML2 = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[30]",))

# Field for Error
CERROR = LocatedComponents(
    phys=PHY.ERRE_R,
    type="ELEM",
    components=(
        "ERREST",
        "NUEST",
        "SIGCAL",
        "TERMRE",
        "TERMR2",
        "TERMNO",
        "TERMN2",
        "TERMSA",
        "TERMS2",
        "TAILLE",
    ),
)

# Field for frequency
CFREQR = LocatedComponents(phys=PHY.FREQ_R, type="ELEM", components=("FREQ",))

# Field for time
MTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))

# Field for material orientation
CCAMASS = LocatedComponents(
    phys=PHY.CAMA_R, type="ELEM", components=("C", "ALPHA", "BETA", "KAPPA", "X", "Y", "Z")
)

# For energy (isoparametric elements)
CENEISO = LocatedComponents(phys=PHY.ENER_R, type="ELEM", components=("TOTALE",))

# For TOU_INI_ELEM
CNTINIF = LocatedComponents(phys=PHY.NEUT_F, type="ELEM", components=("X[30]",))
CNTINIR = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[30]",))
CNTINII = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X[30]",))

# Field for coefficients of norm
CNORMCF = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[30]",))

# Field for flux in thermic (2D - real)
CFLUX2R = LocatedComponents(phys=PHY.FLUX_R, type="ELEM", components=("FLUX", "FLUY"))

# Field for flux in thermic (3D - real)
CFLUX3R = LocatedComponents(phys=PHY.FLUX_R, type="ELEM", components=("FLUX", "FLUY", "FLUZ"))

# Field for flux in thermic (2D - function)
CFLUX2F = LocatedComponents(phys=PHY.FLUX_F, type="ELEM", components=("FLUX", "FLUY"))

# Field for flux in thermic (3D - function)
CFLUX3F = LocatedComponents(phys=PHY.FLUX_F, type="ELEM", components=("FLUX", "FLUY", "FLUZ"))

# Field for time parameters in thermics
CTIMETR = LocatedComponents(
    phys=PHY.INST_R, type="ELEM", components=("INST", "DELTAT", "THETA", "KHI", "R", "RHO")
)

# Field for section (equivalent problem for beams)
CCASECT = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[9]",))

# Field for normal flux in thermic (function)
CFLUXNF = LocatedComponents(phys=PHY.FLUN_F, type="ELEM", components=("FLUN",))

# Field for normal flux in thermic (real)
CFLUXNR = LocatedComponents(phys=PHY.FLUN_R, type="ELEM", components=("FLUN",))

# Field for normal exchange in thermic (coefficient/function)
CHECHPF = LocatedComponents(phys=PHY.COEH_F, type="ELEM", components=("H",))

# Field for normal exchange in thermic (coefficient/real)
CHECHPR = LocatedComponents(phys=PHY.COEH_R, type="ELEM", components=("H",))

# Field for RESI_REFE_RELA COMPONENT NAMES (ALL COMPONENTS, CONSISTENT WITH CRESIREF)
CRESICMP = LocatedComponents(phys=PHY.RESICMP, type="ELEM", components=("X[11]",))

# Field for RESI_REFE_RELA VALUES (ALL COMPONENTS)
CRESIREF = LocatedComponents(
    phys=PHY.RESIREF,
    type="ELEM",
    components=(
        "SIGM",
        "EPSI",
        "FLUXTHER",
        "FLUXHYD[2]",
        "EFFORT",
        "MOMENT",
        "VARI",
        "DEPL",
        "LAGR",
        "PI",
    ),
)

# Field for normal exchange in thermic (temperature/function)
CTEMPEF = LocatedComponents(phys=PHY.TEMP_F, type="ELEM", components=("TEMP",))

# Field for normal exchange in thermic (temperature/real)
CTEMPER = LocatedComponents(phys=PHY.TEMP_R, type="ELEM", components=("TEMP",))

# Field for time parameters in metallurgy
CTIMMTR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST", "DELTA01", "DELTA12"))

# For phases in metallurgy
CPHASES = LocatedComponents(phys=PHY.VAR2_R, type="ELEM", components=("V[9]",))

CTEREFE = LocatedComponents(phys=PHY.TEMP_R, type="ELEM", components=("TEMP",))

# Field for error code
ECODRET = LocatedComponents(phys=PHY.CODE_I, type="ELEM", components=("IRET",))

# Field for values of indicator (AFFE_MODELE)
CINDICR = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[2]",))

# Field for parameters of indicator (AFFE_MODELE)
CCHCKPR = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[2]",))

# Field for parameters of mechanical plate (DKT/DKTG/DST/Q4G/Q4GG)
CSHLDKT = LocatedComponents(
    phys=PHY.CACOQU_R, type="ELEM", components=("EP", "ALPHA", "BETA", "CTOR", "EXCENT", "INERTIE")
)

# Field for parameters of mechanical plate (MITC)
CSHLMITC = LocatedComponents(phys=PHY.CACOQU_R, type="ELEM", components=("EP",))

# Field for parameters of mechanical plate (COQUE_3D)
CSHL3D = LocatedComponents(
    phys=PHY.CACOQU_R,
    type="ELEM",
    components=("EP", "ALPHA", "BETA", "KAPPA", "CTOR", "EXCENT", "INERTIE"),
)

# Field for parameters of mechanical plate (COQUE_AXIS)
CSHLAXI = LocatedComponents(phys=PHY.CACOQU_R, type="ELEM", components=("EP", "KAPPA", "C_METR"))

# Field for parameters of grid (GRILLE_EXCENTRE/GRILLE_MEMBRANE)
CGRID = LocatedComponents(
    phys=PHY.CACOQU_R, type="ELEM", components=("SECT_L", "ALPHA", "BETA", "DIST_N", "CTOR")
)

# Field for parameters of membrane (MEMBRANE)
CMEMB = LocatedComponents(
    phys=PHY.CACOQU_R, type="ELEM", components=("EP", "ALPHA", "BETA", "TENS")
)

# Field for parameters of thermal plate (COQUE_AXIS/COQUE_PLAN)
CSHLTH1 = LocatedComponents(phys=PHY.CACOQU_R, type="ELEM", components=("EP",))

# Field for parameters of thermal plate (COQUE)
CSHLTH2 = LocatedComponents(phys=PHY.CACOQU_R, type="ELEM", components=("EP", "ALPHA", "BETA"))


CNEUTR1 = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X1",))

CBORNPI = LocatedComponents(phys=PHY.PILO_R, type="ELEM", components=("A0", "A1"))

CCACABL = LocatedComponents(phys=PHY.CACABL_R, type="ELEM", components=("SECT", "TENS"))

CCAGNBA = LocatedComponents(phys=PHY.CAGNBA_R, type="ELEM", components=("A1",))

CCAGNP1 = LocatedComponents(
    phys=PHY.CAGNPO_R,
    type="ELEM",
    components=("A1", "IY1", "IZ1", "AY1", "AZ1", "EY1", "EZ1", "JX1", "JG1", "IYR21", "IZR21"),
)

CCAGRPO = LocatedComponents(
    phys=PHY.CAGEPO_R,
    type="ELEM",
    components=(
        "HY1",
        "HZ1",
        "EPY1",
        "EPZ1",
        "HY2",
        "HZ2",
        "EPY2",
        "EPZ2",
        "R1",
        "EP1",
        "R2",
        "EP2",
        "TSEC",
    ),
)

CCARCRI = LocatedComponents(
    phys=PHY.CARCRI,
    type="ELEM",
    components=(
        "ITECREL",
        "MACOMP",
        "RESCREL",
        "THETA",
        "ITEDEC",
        "INTLOC",
        "PERTURB",
        "TOLDEBO",
        "ITEDEBO",
        "RESIRADI",
        "VARIEXT1",
        "THETATHM",
        "POSTITER",
        "LC_EXT[3]",
        "MATRNSYM",
        "ALPHATHM",
        "LC_EXT2[2]",
        "STRAIN",
        "VARIEXT2",
    ),
)

CCDTAU = LocatedComponents(phys=PHY.PILO_R, type="ELEM", components=("A0",))

CCINFDI = LocatedComponents(
    phys=PHY.CINFDI_R,
    type="ELEM",
    components=(
        "REPK",
        "REPM",
        "REPA",
        "SYMK",
        "SYMM",
        "SYMA",
        "DISK",
        "DISM",
        "DISA",
        "ETAK",
        "TYDI",
    ),
)

CCOEFC = LocatedComponents(phys=PHY.IMPE_C, type="ELEM", components=("IMPE",))

CCOEFR = LocatedComponents(phys=PHY.IMPE_R, type="ELEM", components=("IMPE",))

CCOMPOR = LocatedComponents(
    phys=PHY.COMPOR,
    type="ELEM",
    components=(
        "RELCOM",
        "NBVARI",
        "DEFORM",
        "INCELA",
        "C_PLAN",
        "NUME_LC",
        "MULTCOMP",
        "POSTITER",
        "KIT1NAME",
        "KIT2NAME",
        "KIT3NAME",
        "KIT4NAME",
        "KIT1NUME",
        "KIT2NUME",
        "KIT3NUME",
        "KIT4NUME",
        "KIT1NVAR",
        "KIT2NVAR",
        "KIT3NVAR",
        "KIT4NVAR",
        "DEFO_LDC",
        "RIGIGEOM",
        "REGUVISC",
        "MGISADDR",
        "POSTINCR",
    ),
)

CCOMPOT = LocatedComponents(
    phys=PHY.COMPOR, type="ELEM", components=("RELCOM", "NBVARI", "DEFORM", "INCELA", "C_PLAN")
)

CCMBHHO = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[8]",))

CMPHHOR = LocatedComponents(phys=PHY.N480_R, type="ELEM", components=("X[240]",))

CCOMPO2 = LocatedComponents(phys=PHY.COMPOR, type="ELEM", components=("NBVARI",))

CMLCOMP = LocatedComponents(phys=PHY.MULTCOMP, type="ELEM", components=("SD_COMP",))

CCONFR = LocatedComponents(phys=PHY.N120_R, type="ELEM", components=("X[60]",))

CCONPT = LocatedComponents(phys=PHY.N120_R, type="ELEM", components=("X[34]",))

CCONSTR = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X1",))

CDECENT = LocatedComponents(phys=PHY.NEUT_K24, type="ELEM", components=("Z1",))

CFER1_R = LocatedComponents(
    phys=PHY.FER1_R,
    type="ELEM",
    components=(
        "TYPCOMB",
        "CODIF",
        "METH2D",
        "THITER",
        "EPITER",
        "APHITER",
        "COND109",
        "TYPSTRU",
        "FERRSYME",
        "SLSYME",
        "FERRCOMP",
        "EPUCISA",
        "FERRMIN",
        "RHOLMIN",
        "RHOTMIN",
        "COMPRESS",
        "CEQUI",
        "ENROBI",
        "ENROBS",
        "ENROBYI",
        "ENROBYS",
        "ENROBZI",
        "ENROBZS",
        "SIGS",
        "SIGCI",
        "SIGCS",
        "SIGCYI",
        "SIGCYS",
        "SIGCZI",
        "SIGCZS",
        "ALPHACC",
        "GAMMAS",
        "GAMMAC",
        "FACIER",
        "EYS",
        "TYPDIAG",
        "FBETON",
        "CLACIER",
        "UC",
        "UM",
        "RHOACIER",
        "AREINF",
        "ASHEAR",
        "ASTIRR",
        "RHOCRIT",
        "DATCRIT",
        "LCRIT",
        "WMAXI",
        "WMAXS",
        "WMAXYI",
        "WMAXYS",
        "WMAXZI",
        "WMAXZS",
        "SIGELSQP",
        "KT",
        "PHIXI",
        "PHIXS",
        "PHIYI",
        "PHIYS",
        "PHIZI",
        "PHIZS",
        "ASMINFIS",
        "ECHELLE",
    ),
)

CFER2_R = LocatedComponents(
    phys=PHY.FER2_R,
    type="ELEM",
    components=(
        "DNSXI",
        "DNSXS",
        "DNSYI",
        "DNSYS",
        "DNSXT",
        "DNSYT",
        "DNSVOL",
        "CONSTRUC",
        "ASMINXI",
        "ASMINXS",
        "ASMINYI",
        "ASMINYS",
    ),
)

CVFER0_R = LocatedComponents(
    phys=PHY.FER2_R, type="ELEM", components=("DNSXI", "DNSXS", "DNSYI", "DNSYS")
)


CVFER1_R = LocatedComponents(
    phys=PHY.VFER1_R,
    type="ELEM",
    components=(
        "TYPCOMB",
        "CODIF",
        "THITER",
        "TYPSTRU",
        "ENROBI",
        "ENROBS",
        "SIGS",
        "SIGCI",
        "SIGCS",
        "ALPHACC",
        "GAMMAS",
        "GAMMAC",
        "FACIER",
        "EYS",
        "TYPDIAG",
        "FBETON",
        "CLACIER",
        "UC",
        "UM",
        "CEQUI",
    ),
)

CVFER2_R = LocatedComponents(
    phys=PHY.VFER2_R,
    type="ELEM",
    components=(
        "MARGE",
        "THETA",
        "TAU_UTIL",
        "C0C",
        "C0CRD",
        "DNSINF",
        "DNSSUP",
        "N0",
        "M0",
        "NED",
        "MED",
        "N_DIAG",
        "M_DIAG",
        "TYPCMB",
        "TYPCO",
        "CEQUI",
        "ENROBI",
        "ENROBS",
        "SIGS",
        "SIGCI",
        "SIGCS",
        "ALPHACC",
        "GAMMAS",
        "GAMMAC",
        "FACIER",
        "EYS",
        "TYPDIAG",
        "FBETON",
        "CLACIER",
        "UC",
        "UM",
        "HT",
        "BW",
    ),
)

CFER3_R = LocatedComponents(
    phys=PHY.FER2_R,
    type="ELEM",
    components=("AYI", "AYS", "AZI", "AZS", "AST", "ATOT", "DNSVOL", "CONSTRUC"),
)

CFISSR = LocatedComponents(
    phys=PHY.FISS_R, type="ELEM", components=("XA", "YA", "XTAN", "YTAN", "XNORM", "YNORM")
)

CFLAPLA = LocatedComponents(phys=PHY.FLAPLA, type="ELEM", components=("NOMAIL", "NOGEOM"))

CFORCEC = LocatedComponents(
    phys=PHY.FORC_C, type="ELEM", components=("FX", "FY", "FZ", "MX", "MY", "MZ", "REP")
)

CFTRC = LocatedComponents(phys=PHY.ADRSJEVN, type="ELEM", components=("I[2]",))

CHARMON = LocatedComponents(phys=PHY.HARMON, type="ELEM", components=("NH",))


CIMPEDC = LocatedComponents(phys=PHY.IMPE_C, type="ELEM", components=("IMPE",))

CINSTPR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))

CITERAT = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X1",))

CLISTMA = LocatedComponents(phys=PHY.NEUT_K16, type="ELEM", components=("Z[2]",))

CMASDIA = LocatedComponents(phys=PHY.POSI, type="ELEM", components=("POS",))

CMATERC = LocatedComponents(phys=PHY.ADRSJEVE, type="ELEM", components=("I1",))

CNUMMOD = LocatedComponents(phys=PHY.NUMMOD, type="ELEM", components=("NUM",))

COMEG2R = LocatedComponents(phys=PHY.OME2_R, type="ELEM", components=("OMEG2",))

CONDPLA = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z[2]",))

CONDPLR = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[10]",))

CPAIRR = LocatedComponents(phys=PHY.N120_R, type="ELEM", components=("X[27]",))

CPESANR = LocatedComponents(phys=PHY.PESA_R, type="ELEM", components=("G", "AG", "BG", "CG"))

CPHASES = LocatedComponents(phys=PHY.VAR2_R, type="ELEM", components=("V[9]",))

CRAYONF = LocatedComponents(phys=PHY.RAYO_F, type="ELEM", components=("SIGMA", "EPSIL", "TPINF"))

CRAYONR = LocatedComponents(phys=PHY.RAYO_R, type="ELEM", components=("SIGMA", "EPSIL", "TPINF"))

CREFERI = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X[12]",))

CREFERK = LocatedComponents(phys=PHY.NEUT_K24, type="ELEM", components=("Z[19]",))

CROTATR = LocatedComponents(
    phys=PHY.ROTA_R, type="ELEM", components=("OME", "AR", "BR", "CR", "X", "Y", "Z")
)

CSOURCF = LocatedComponents(phys=PHY.SOUR_F, type="ELEM", components=("SOUR",))

CSOURCR = LocatedComponents(phys=PHY.SOUR_R, type="ELEM", components=("SOUR",))

CSOUSOP = LocatedComponents(phys=PHY.NEUT_K24, type="ELEM", components=("Z1",))

CSUROPT = LocatedComponents(phys=PHY.NEUT_K24, type="ELEM", components=("Z1",))


CTYPEPI = LocatedComponents(phys=PHY.PILO_K, type="ELEM", components=("TYPE",))

CT_EXTR = LocatedComponents(
    phys=PHY.TEMP_R, type="ELEM", components=("TEMP", "TEMP_INF", "TEMP_SUP")
)

CVENTCX = LocatedComponents(phys=PHY.VENTCX_F, type="ELEM", components=("FCXP",))

E120NEUT = LocatedComponents(phys=PHY.N816_R, type="ELEM", components=("X[120]",))

E10NEUTI = LocatedComponents(phys=PHY.N120_I, type="ELEM", components=("X[10]",))

E128NEUI = LocatedComponents(phys=PHY.N512_I, type="ELEM", components=("X[128]",))

E12NEUTR = LocatedComponents(phys=PHY.N132_R, type="ELEM", components=("X[12]",))

E132NEUR = LocatedComponents(phys=PHY.N132_R, type="ELEM", components=("X[132]",))

E144NEUI = LocatedComponents(phys=PHY.N1280I, type="ELEM", components=("X[144]",))

E14NEUTR = LocatedComponents(phys=PHY.N816_R, type="ELEM", components=("X[14]",))

E15NEUTR = LocatedComponents(phys=PHY.N480_R, type="ELEM", components=("X[15]",))

E162NEUR = LocatedComponents(phys=PHY.N2448R, type="ELEM", components=("X[162]",))

E170NEUT = LocatedComponents(phys=PHY.N1360R, type="ELEM", components=("X[170]",))

E200NEUT = LocatedComponents(phys=PHY.N1360R, type="ELEM", components=("X[200]",))

E18NEUI = LocatedComponents(phys=PHY.N120_I, type="ELEM", components=("X[18]",))

E198NEUT = LocatedComponents(phys=PHY.N792_R, type="ELEM", components=("X[198]",))

E1NEUK8 = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z1",))

E1NEUTI = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X1",))

E1NEUTR = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X1",))

E220NEUR = LocatedComponents(phys=PHY.N480_R, type="ELEM", components=("X[220]",))

E22NEUTR = LocatedComponents(phys=PHY.N792_R, type="ELEM", components=("X[22]",))

E24NEUTR = LocatedComponents(phys=PHY.N132_R, type="ELEM", components=("X[24]",))

E28NEUTR = LocatedComponents(phys=PHY.N2448R, type="ELEM", components=("X[28]",))

E2NEUTI = LocatedComponents(phys=PHY.N512_I, type="ELEM", components=("X[2]",))

E20NEUTI = LocatedComponents(phys=PHY.N120_I, type="ELEM", components=("X[20]",))

E2NEUTR = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[2]",))

E306NEUT = LocatedComponents(phys=PHY.N2448R, type="ELEM", components=("X[306]",))

E360NEUT = LocatedComponents(phys=PHY.N2448R, type="ELEM", components=("X[360]",))

E320NEUI = LocatedComponents(phys=PHY.N1280I, type="ELEM", components=("X[320]",))

E32NEUTI = LocatedComponents(phys=PHY.N512_I, type="ELEM", components=("X[32]",))

E32NEUTR = LocatedComponents(phys=PHY.N2448R, type="ELEM", components=("X[32]",))

E35NEUTR = LocatedComponents(phys=PHY.N1360R, type="ELEM", components=("X[35]",))

E36NEUI = LocatedComponents(phys=PHY.N1280I, type="ELEM", components=("X[36]",))

E3NEUTI = LocatedComponents(phys=PHY.N120_I, type="ELEM", components=("X[3]",))

E40NEUTR = LocatedComponents(phys=PHY.N1360R, type="ELEM", components=("X[40]",))

E4NEUTR = LocatedComponents(phys=PHY.N792_R, type="ELEM", components=("X[4]",))

E512NEUI = LocatedComponents(phys=PHY.N1280I, type="ELEM", components=("X[512]",))

E54NEUTI = LocatedComponents(phys=PHY.N720_I, type="ELEM", components=("X[54]",))

E54NEUTR = LocatedComponents(phys=PHY.N816_R, type="ELEM", components=("X[54]",))

E55NEUTR = LocatedComponents(phys=PHY.N480_R, type="ELEM", components=("X[55]",))

E60NEUTR = LocatedComponents(phys=PHY.N480_R, type="ELEM", components=("X[60]",))

E5NEUTR = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[5]",))

E6NEUTR = LocatedComponents(phys=PHY.N132_R, type="ELEM", components=("X[6]",))

E72NEUI = LocatedComponents(phys=PHY.N1280I, type="ELEM", components=("X[72]",))

E88NEUTR = LocatedComponents(phys=PHY.N792_R, type="ELEM", components=("X[88]",))

E8NEUI = LocatedComponents(phys=PHY.N960_I, type="ELEM", components=("X[8]",))

E90NEUTI = LocatedComponents(phys=PHY.N720_I, type="ELEM", components=("X[90]",))

E9NEUTI = LocatedComponents(phys=PHY.N720_I, type="ELEM", components=("X[9]",))

ECAFIEL = LocatedComponents(
    phys=PHY.CAFI_R, type="ELEM", components=("YG", "ZG", "AIRE", "YP", "ZP", "GX", "NUMGR")
)

ECARAGE = LocatedComponents(
    phys=PHY.MASS_R,
    type="ELEM",
    components=(
        "M",
        "CDGX",
        "CDGY",
        "CDGZ",
        "IXX",
        "IYY",
        "IZZ",
        "IXY",
        "IXZ",
        "IYZ",
        "IXR2",
        "IYR2",
    ),
)

ECHALIM = LocatedComponents(phys=PHY.CHLI_R, type="ELEM", components=("CHLI[3]",))


ECOOR1R = LocatedComponents(phys=PHY.N120_R, type="ELEM", components=("X[81]",))

ECOURAN = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X1",))

EDCEL_I = LocatedComponents(phys=PHY.DCEL_I, type="ELEM", components=("NPG_DYN", "NCMP_DYN"))

EEMATE_R = LocatedComponents(
    phys=PHY.MATE_R,
    type="ELEM",
    components=("X", "Y", "Z", "E", "NU", "RHO", "ALPHA", "LAMBDA", "RHO_CP"),
)

EENECNO = LocatedComponents(
    phys=PHY.ENER_R, type="ELEM", components=("TOTALE", "DX", "DY", "DZ", "DRX", "DRY", "DRZ")
)

EENEDNO = LocatedComponents(
    phys=PHY.ENER_R, type="ELEM", components=("TOTALE", "TRAC_COM", "TORSION", "FLEX_Y", "FLEX_Z")
)

EERREURT = LocatedComponents(
    phys=PHY.ERRE_R,
    type="ELEM",
    components=(
        "ERTABS",
        "ERTREL",
        "TERMNO",
        "TERMVO",
        "TERMV2",
        "TERMV1",
        "TERMSA",
        "TERMS2",
        "TERMS1",
        "TERMFL",
        "TERMF2",
        "TERMF1",
        "TERMEC",
        "TERME2",
        "TERME1",
    ),
)

ETHETA = LocatedComponents(
    phys=PHY.THET_R,
    type="ELNO",
    components=("MODULE", "DIR_X", "DIR_Y", "DIR_Z", "ABSC_CUR", "LONG"),
)


EHECHPR = LocatedComponents(phys=PHY.COEH_R, type="ELEM", components=("H",))

EMASSINE = LocatedComponents(
    phys=PHY.MASS_R,
    type="ELEM",
    components=("M", "CDGX", "CDGY", "CDGZ", "IXX", "IYY", "IZZ", "IXY", "IXZ", "IYZ"),
)

EMNEUT_I = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X1",))

ENONLIN = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X1",))

ENORME = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X1",))

EPJSIGM = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELEM",
    components=(
        "SIG_NX",
        "SIG_NY",
        "SIG_NZ",
        "SIG_N",
        "SIG_TX",
        "SIG_TY",
        "SIG_TZ",
        "SIG_T1X",
        "SIG_T1Y",
        "SIG_T1Z",
        "SIG_T1",
        "SIG_T2X",
        "SIG_T2Y",
        "SIG_T2Z",
        "SIG_T2",
        "SIG_TN",
    ),
)

EREFE1K = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z1",))

ERICTRA = LocatedComponents(
    phys=PHY.RICE_TRA, type="ELEM", components=("TRIAX", "RSR0", "VOLU", "NUMEMA", "DEPSEQ")
)

ESINGUL = LocatedComponents(phys=PHY.SING_R, type="ELEM", components=("DEGRE", "RAPPORT", "TAILLE"))

ET_EXTR = LocatedComponents(phys=PHY.TEMP_R, type="ELEM", components=("TEMP",))

EVOISIN = LocatedComponents(phys=PHY.VOISIN, type="ELEM", components=("V0", "V[6]", "T0", "T[6]"))

EWEIBUL = LocatedComponents(phys=PHY.WEIBULL, type="ELEM", components=("DSIGWB",))

FISCO_I = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X[2]",))

I1NEUT_I = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X1",))

I3NEUT_I = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X[3]",))

MALPHAR = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X1",))

MDDLIMC = LocatedComponents(phys=PHY.DDLI_C, type="ELEM", components=("C",))

MDDLIMF = LocatedComponents(phys=PHY.DDLI_F, type="ELEM", components=("C",))

MDDLIMR = LocatedComponents(phys=PHY.DDLI_R, type="ELEM", components=("C",))

MDDLMUR = LocatedComponents(phys=PHY.DDLM_R, type="ELEM", components=("A1",))

N40NEUI = LocatedComponents(phys=PHY.N120_I, type="ELEM", components=("X[40]",))

N80NEUI = LocatedComponents(phys=PHY.N120_I, type="ELEM", components=("X[80]",))

NFAMILK = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z1",))

NINFORR = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[26]",))

# Field for load VITE_FACE (complex)
EVITEFC = LocatedComponents(
    phys=PHY.VFAC_C, type="ELEM", components=("VITE", "INDC", "DIRX", "DIRY", "DIRZ")
)

# Field for load VITE_FACE (function)
EVITEFF = LocatedComponents(
    phys=PHY.VFAC_F, type="ELEM", components=("VITE", "INDC", "DIRX", "DIRY", "DIRZ")
)

# Field for load VITE_FACE (real)
EVITEFR = LocatedComponents(
    phys=PHY.VFAC_R, type="ELEM", components=("VITE", "INDC", "DIRX", "DIRY", "DIRZ")
)

# ----------------------------------------------------------------------------------------------
# Located components - ELNO - Field on nodes by element
# ----------------------------------------------------------------------------------------------
# For TOU_INI_ELNO
ENTINIF = LocatedComponents(phys=PHY.NEUT_F, type="ELNO", components=("X[30]",))
ENTINIR = LocatedComponents(phys=PHY.NEUT_R, type="ELNO", components=("X[30]",))
ENTINII = LocatedComponents(phys=PHY.NEUT_I, type="ELNO", components=("X[30]",))


# For stresses (Real-3D)
ESIG3DR = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
)

# For stresses (Complex-3D)
ESIG3DC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
)

# For stresses (Real-2D)
ESIG2DR = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY")
)

# For stresses (Complex-2D)
ESIG2DC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY")
)

# For geometry (2D)
EGEOM2D = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))

# For geometry (3D)
EGEOM3D = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))

# For pressure (Real-3D)
EPRE3DR = LocatedComponents(phys=PHY.PRES_R, type="ELNO", components=("PRES",))

# For pressure (Real-3D)
EPRE2DR = LocatedComponents(phys=PHY.PRES_R, type="ELNO", components=("PRES", "CISA"))

# For equivalent stresses (Real-3D)
ECOEQNO = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELNO",
    components=(
        "VMIS",
        "TRESCA",
        "PRIN_[3]",
        "VMIS_SG",
        "VECT_1_X",
        "VECT_1_Y",
        "VECT_1_Z",
        "VECT_2_X",
        "VECT_2_Y",
        "VECT_2_Z",
        "VECT_3_X",
        "VECT_3_Y",
        "VECT_3_Z",
        "TRSIG",
        "TRIAX",
    ),
)

# For strains (Real-3D)
EEPS3DR = LocatedComponents(
    phys=PHY.EPSI_R, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)

# For strains (Complex-3D)
EEPS3DC = LocatedComponents(
    phys=PHY.EPSI_C, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)

# For strains (Real-2D)
EEPS2DR = LocatedComponents(
    phys=PHY.EPSI_R, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)

# For strains (Complex-2D)
EEPS2DC = LocatedComponents(
    phys=PHY.EPSI_C, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)

# For equivalent strains (Real-3D)
EDFEQNO = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELNO",
    components=(
        "INVA_2",
        "PRIN_[3]",
        "INVA_2SG",
        "VECT_1_X",
        "VECT_1_Y",
        "VECT_1_Z",
        "VECT_2_X",
        "VECT_2_Y",
        "VECT_2_Z",
        "VECT_3_X",
        "VECT_3_Y",
        "VECT_3_Z",
    ),
)

# For external state variables strains (3D)
NVARC3D = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELNO",
    components=("EPTHER_L", "EPTHER_T", "EPTHER_N", "EPSECH", "EPHYDR", "EPPTOT"),
)

# For external state variables strains (1D)
NVARC1D = LocatedComponents(phys=PHY.EPSI_R, type="ELNO", components=("EPTHER_L",))

# For nodal forces (Real-2D)
NFOR2DR = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY"))

# For nodal forces (Real-3D)
NFOR3DR = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY", "FZ"))

# For energy (isoparametric elements)
NENEISO = LocatedComponents(phys=PHY.ENER_R, type="ELNO", components=("TOTALE",))

# For error
NERROR = LocatedComponents(
    phys=PHY.ERRE_R,
    type="ELNO",
    components=(
        "ERREST",
        "NUEST",
        "SIGCAL",
        "TERMRE",
        "TERMR2",
        "TERMNO",
        "TERMN2",
        "TERMSA",
        "TERMS2",
        "TAILLE",
    ),
)

# For speed in thermic (2D)
NVITE2R = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY"))

# For speed in thermic (3D)
NVITE3R = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))

# Field for flux in thermic (2D - real)
NFLUX2R = LocatedComponents(phys=PHY.FLUX_R, type="ELNO", components=("FLUX", "FLUY"))

# Field for flux in thermic (3D - real)
NFLUX3R = LocatedComponents(phys=PHY.FLUX_R, type="ELNO", components=("FLUX", "FLUY", "FLUZ"))

# Field for gradient in thermic (2D - real)
NGRAT2R = LocatedComponents(phys=PHY.GRAT_R, type="ELNO", components=("GRAT_X", "GRAT_Y"))

# Field for gradient in thermic (3D - real)
NGRAT3R = LocatedComponents(phys=PHY.GRAT_R, type="ELNO", components=("GRAT_X", "GRAT_Y", "GRAT_Z"))

# For body source in thermic
NSOURCR = LocatedComponents(phys=PHY.SOUR_R, type="ELNO", components=("SOUR",))

# For phases in metallurgy
EPHASES = LocatedComponents(phys=PHY.VARI_R, type="ELNO", components=("VARI",))

# For generalized forces - Beams (Real)
NEFGEBR = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("N", "VY", "VZ", "MT", "MFY", "MFZ")
)

# For generalized forces - Beams (Complex)
NEFGEBC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("N", "VY", "VZ", "MT", "MFY", "MFZ")
)

# For MIN-MAX SP
NMINMAX = LocatedComponents(
    phys=PHY.SPMX_R, type="ELNO", components=("VAL", "NUCOU", "NUSECT", "NUFIBR", "POSIC", "POSIS")
)

DDL_NOZ1 = LocatedComponents(phys=PHY.SIZZ_R, type="ELNO", components=("SIZZ",))

E3NEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELNO", components=("X[3]",))

ECONTPO = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("SN", "SVY", "SVZ", "SMT", "SMFY", "SMFZ")
)

ECOPILN = LocatedComponents(phys=PHY.PILO_R, type="ELNO", components=("A0", "A[3]", "ETA"))

EDEFGNC = LocatedComponents(
    phys=PHY.EPSI_C,
    type="ELNO",
    components=("EXX", "EYY", "EXY", "KXX", "KYY", "KXY", "GAX", "GAY"),
)

EDERANO = LocatedComponents(
    phys=PHY.DERA_R,
    type="ELNO",
    components=(
        "DCHA_V",
        "DCHA_T",
        "IND_DCHA",
        "VAL_DCHA",
        "X11",
        "X22",
        "X33",
        "X12",
        "X13",
        "X23",
        "RADI_V",
        "ERR_RADI",
    ),
)

EDIFFNO = LocatedComponents(phys=PHY.DIFF_R, type="ELNO", components=("DIFT", "DIFL", "DIFV"))

EDOMGNO = LocatedComponents(phys=PHY.DOMA_R, type="ELNO", components=("DOMA",))

EDURTNO = LocatedComponents(phys=PHY.DURT_R, type="ELNO", components=("HV",))

EEFGENOQ = LocatedComponents(phys=PHY.SIEF_R, type="ELNO", components=("MT", "MFY", "MFZ", "MEQ"))

EEINST_R = LocatedComponents(phys=PHY.INST_R, type="ELNO", components=("INST",))


EERRENOT = LocatedComponents(
    phys=PHY.ERRE_R,
    type="ELNO",
    components=(
        "ERTABS",
        "ERTREL",
        "TERMNO",
        "TERMVO",
        "TERMV2",
        "TERMV1",
        "TERMSA",
        "TERMS2",
        "TERMS1",
        "TERMFL",
        "TERMF2",
        "TERMF1",
        "TERMEC",
        "TERME2",
        "TERME1",
    ),
)

EHYDRNO = LocatedComponents(phys=PHY.HYDR_R, type="ELNO", components=("HYDR",))

EHYGRNO = LocatedComponents(phys=PHY.HYGR_R, type="ELNO", components=("HUMR", "PCAP"))

ENINST_R = LocatedComponents(phys=PHY.INST_R, type="ELNO", components=("INST",))

EPRACNO = LocatedComponents(phys=PHY.PRAC_R, type="ELNO", components=("PRES_R", "PRES_I", "DB"))

EPRMENO = LocatedComponents(phys=PHY.PRME_R, type="ELNO", components=("DB",))

ESIMXNC = LocatedComponents(phys=PHY.SIEFMX_C, type="ELNO", components=("SIXXMIN", "SIXXMAX"))

ESIMXNO = LocatedComponents(phys=PHY.SIEFMX_R, type="ELNO", components=("SIXXMIN", "SIXXMAX"))

ESINGNO = LocatedComponents(phys=PHY.SING_R, type="ELNO", components=("DEGRE", "RAPPORT", "TAILLE"))

ETRIANO = LocatedComponents(
    phys=PHY.ENDO_R, type="ELNO", components=("TRIAX", "SI_ENDO", "COENDO", "DOM_LEM")
)

FISNO_I = LocatedComponents(phys=PHY.NEUT_I, type="ELNO", components=("X1",))

N1NEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELNO", components=("X1",))

N2NEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELNO", components=("X[2]",))

N3NEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELNO", components=("X[3]",))

N5NEUTR = LocatedComponents(phys=PHY.NEUT_R, type="ELNO", components=("X[5]",))

N5NEUTI = LocatedComponents(phys=PHY.N120_I, type="ELNO", components=("X[5]",))

N6NEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELNO", components=("X[6]",))

N9NEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELNO", components=("X[9]",))

NTEMPER = LocatedComponents(
    phys=PHY.TEMP_R, type="ELNO", components=("TEMP_MIL", "TEMP_INF", "TEMP_SUP")
)

ZVARINO = LocatedComponents(phys=PHY.VARI_R, type="ELNO", components=("VARI",))

ZVARCNO = LocatedComponents(phys=PHY.VARI_R, type="ELNO", components=("VARI",))

# ----------------------------------------------------------------------------------------------
# Located components - ELGA - Field on integration points
# ----------------------------------------------------------------------------------------------
# Coordinates/weight of Gauss points (in 2D)
EGGAU2D = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "W")
)

# Coordinates/weight of Gauss points (in 3D)
EGGAU3D = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)

# For TOU_INI_ELGA
EGTINIF = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))
EGTINIR = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))
EGTINII = LocatedComponents(phys=PHY.NEUT_I, type="ELGA", location="RIGI", components=("X[30]",))

EGTINIV = LocatedComponents(
    phys=PHY.VALO_R, type="ELGA", location="RIGI", components=("VALEUR", "GRAD_X", "GRAD_Y")
)

# Geometry at Gauss points (in 2D)
EGGEO2D = LocatedComponents(phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y"))

# Geometry at Gauss points (in 3D)
EGGEO3D = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z")
)

# Displacements at Gauss points (in 2D)
EGDEP2D = LocatedComponents(phys=PHY.DEPL_R, type="ELGA", location="RIGI", components=("DX", "DY"))

# Displacements at Gauss points (in 3D)
EGDEP3D = LocatedComponents(
    phys=PHY.DEPL_R, type="ELGA", location="RIGI", components=("DX", "DY", "DZ")
)

# For stresses (Real-3D)
EGIG3DR = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"),
)

# For stresses (Complex-3D)
EGIG3DC = LocatedComponents(
    phys=PHY.SIEF_C,
    type="ELGA",
    location="RIGI",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"),
)

# For stresses (Real-2D)
EGIG2DR = LocatedComponents(
    phys=PHY.SIEF_R, type="ELGA", location="RIGI", components=("SIXX", "SIYY", "SIZZ", "SIXY")
)

# For stresses (Complex-2D)
EGIG2DC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELGA", location="RIGI", components=("SIXX", "SIYY", "SIZZ", "SIXY")
)

# For equivalent stresses (Real)
ECOEQPG = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=(
        "VMIS",
        "TRESCA",
        "PRIN_[3]",
        "VMIS_SG",
        "VECT_1_X",
        "VECT_1_Y",
        "VECT_1_Z",
        "VECT_2_X",
        "VECT_2_Y",
        "VECT_2_Z",
        "VECT_3_X",
        "VECT_3_Y",
        "VECT_3_Z",
        "TRSIG",
        "TRIAX",
    ),
)

# For internal state variables
ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))

# For external state variables
ZVARCPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="MATER", components=("VARI",))

# For strains (Real-3D)
EGPS3DR = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"),
)

# For strains (Complex-3D)
EGPS3DC = LocatedComponents(
    phys=PHY.EPSI_C,
    type="ELGA",
    location="RIGI",
    components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"),
)

# For strains (Real-2D)
EGPS2DR = LocatedComponents(
    phys=PHY.EPSI_R, type="ELGA", location="RIGI", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)

# For strains (Complex-2D)
EGPS2DC = LocatedComponents(
    phys=PHY.EPSI_C, type="ELGA", location="RIGI", components=("EPXX", "EPYY", "EPZZ", "EPXY")
)

# For equivalent strains (Real)
EDFEQPG = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=(
        "INVA_2",
        "PRIN_[3]",
        "INVA_2SG",
        "VECT_1_X",
        "VECT_1_Y",
        "VECT_1_Z",
        "VECT_2_X",
        "VECT_2_Y",
        "VECT_2_Z",
        "VECT_3_X",
        "VECT_3_Y",
        "VECT_3_Z",
    ),
)

# For external state variables strains (3D)
EGVARC3D = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EPTHER_L", "EPTHER_T", "EPTHER_N", "EPSECH", "EPHYDR", "EPPTOT"),
)

# For external state variables strains (1D)
EGVARC1D = LocatedComponents(
    phys=PHY.EPSI_R, type="ELGA", location="RIGI", components=("EPTHER_L",)
)

# For fatigue (3D)
EGFC3DR = LocatedComponents(
    phys=PHY.FACY_R,
    type="ELGA",
    location="RIGI",
    components=(
        "DTAUM1",
        "VNM1X",
        "VNM1Y",
        "VNM1Z",
        "SINMAX1",
        "SINMOY1",
        "EPNMAX1",
        "EPNMOY1",
        "SIGEQ1",
        "NBRUP1",
        "ENDO1",
        "DTAUM2",
        "VNM2X",
        "VNM2Y",
        "VNM2Z",
        "SINMAX2",
        "SINMOY2",
        "EPNMAX2",
        "EPNMOY2",
        "SIGEQ2",
        "NBRUP2",
        "ENDO2",
        "VMIS",
        "TRESCA",
    ),
)

# Field for time
EGINST_R = LocatedComponents(phys=PHY.INST_R, type="ELGA", location="RIGI", components=("INST",))

# For energy (isoparametric elements)
EENEISO = LocatedComponents(phys=PHY.ENER_R, type="ELGA", location="RIGI", components=("TOTALE",))

# For nodal forces (Real-2D)
EFOR2DR = LocatedComponents(phys=PHY.FORC_R, type="ELGA", location="RIGI", components=("FX", "FY"))

# For nodal forces (Real-3D)
EFOR3DR = LocatedComponents(
    phys=PHY.FORC_R, type="ELGA", location="RIGI", components=("FX", "FY", "FZ")
)

# For nodal forces (Real-3D)
ENFOR3R = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY", "FZ"))


# For body source in thermic
ESOURCR = LocatedComponents(phys=PHY.SOUR_R, type="ELGA", location="RIGI", components=("SOUR",))

# For continuation method (PILOTAGE)
ECOPILO = LocatedComponents(
    phys=PHY.PILO_R, type="ELGA", location="RIGI", components=("A0", "A[3]", "ETA")
)

# Field for diffusion in drying (real)
EDIFFUR = LocatedComponents(
    phys=PHY.DIFF_R, type="ELGA", location="RIGI", components=("DIFT", "DIFL", "DIFV")
)

# Field for flux in thermic (2D - real)
EFLUX2R = LocatedComponents(
    phys=PHY.FLUX_R, type="ELGA", location="RIGI", components=("FLUX", "FLUY")
)

# Field for flux in thermic (3D - real)
EFLUX3R = LocatedComponents(
    phys=PHY.FLUX_R, type="ELGA", location="RIGI", components=("FLUX", "FLUY", "FLUZ")
)

# Field for gradient in thermic (2D - real)
EGRAT2R = LocatedComponents(
    phys=PHY.GRAT_R, type="ELGA", location="RIGI", components=("GRAT_X", "GRAT_Y")
)

# Field for gradient in thermic (3D - real)
EGRAT3R = LocatedComponents(
    phys=PHY.GRAT_R, type="ELGA", location="RIGI", components=("GRAT_X", "GRAT_Y", "GRAT_Z")
)

# Field for hydratation
EHYDRR = LocatedComponents(phys=PHY.HYDR_R, type="ELGA", location="MASS", components=("HYDR",))

# Field for hygrometrie in drying (real)
EHYGROR = LocatedComponents(
    phys=PHY.HYGR_R, type="ELGA", location="RIGI", components=("HUMR", "PCAP")
)

# For generalized forces - Plates (Real)
EGFG3DR = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=("NXX", "NYY", "NXY", "MXX", "MYY", "MXY", "QX", "QY"),
)

# For generalized forces - Plates (Complex)
EGFG3DC = LocatedComponents(
    phys=PHY.SIEF_C,
    type="ELGA",
    location="RIGI",
    components=("NXX", "NYY", "NXY", "MXX", "MYY", "MXY", "QX", "QY"),
)

# For generalized forces - Beams (Real)
EEFGEBR = LocatedComponents(
    phys=PHY.SIEF_R, type="ELGA", location="RIGI", components=("N", "VY", "VZ", "MT", "MFY", "MFZ")
)

# For generalized forces - Beams (Complex)
EEFGEBC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELGA", location="RIGI", components=("N", "VY", "VZ", "MT", "MFY", "MFZ")
)

# For MIN-MAX SP
EGMINMAX = LocatedComponents(
    phys=PHY.SPMX_R,
    type="ELGA",
    location="RIGI",
    components=("VAL", "NUCOU", "NUSECT", "NUFIBR", "POSIC", "POSIS"),
)

EGMATE_R = LocatedComponents(
    phys=PHY.MATE_R,
    type="ELGA",
    location="MTGA",
    components=("X", "Y", "Z", "E", "NU", "RHO", "ALPHA", "LAMBDA", "RHO_CP"),
)

ENMATE_R = LocatedComponents(
    phys=PHY.MATE_R,
    type="ELNO",
    components=("X", "Y", "Z", "E", "NU", "RHO", "ALPHA", "LAMBDA", "RHO_CP"),
)


EDERAPG = LocatedComponents(
    phys=PHY.DERA_R,
    type="ELGA",
    location="RIGI",
    components=(
        "DCHA_V",
        "DCHA_T",
        "IND_DCHA",
        "VAL_DCHA",
        "X11",
        "X22",
        "X33",
        "X12",
        "X13",
        "X23",
        "RADI_V",
        "ERR_RADI",
    ),
)

EDOMGGA = LocatedComponents(phys=PHY.DOMA_R, type="ELGA", location="RIGI", components=("DOMA",))


EGGEMA_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="MATER", components=("X", "Y", "Z")
)

EGINDLO = LocatedComponents(
    phys=PHY.INDL_R, type="ELGA", location="RIGI", components=("INDICE", "DIR[4]")
)

EGNEUT1R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X1",))

EHYDRMA = LocatedComponents(phys=PHY.HYDR_R, type="ELGA", location="MATER", components=("HYDR",))

EIMPEDR = LocatedComponents(phys=PHY.IMPE_R, type="ELGA", location="RIGI", components=("IMPE",))

EONDEPR = LocatedComponents(phys=PHY.ONDE_R, type="ELGA", location="RIGI", components=("PRES",))

EPDILPG = LocatedComponents(phys=PHY.PDIL_R, type="ELGA", location="RIGI", components=("A1_LC2",))

EPRESGA = LocatedComponents(phys=PHY.PRES_R, type="ELGA", location="RIGI", components=("PRES",))

ETEMPMA = LocatedComponents(phys=PHY.TEMP_R, type="ELGA", location="MATER", components=("TEMP",))

ETEMPNO = LocatedComponents(phys=PHY.TEMP_R, type="ELNO", components=("TEMP", "SECH"))

ETEMPPG = LocatedComponents(phys=PHY.TEMP_R, type="ELGA", location="RIGI", components=("TEMP",))

ESECHPG = LocatedComponents(phys=PHY.TEMP_R, type="ELGA", location="RIGI", components=("SECH",))

ETRIAPG = LocatedComponents(
    phys=PHY.ENDO_R,
    type="ELGA",
    location="RIGI",
    components=("TRIAX", "SI_ENDO", "COENDO", "DOM_LEM"),
)

EVARC_R = LocatedComponents(
    phys=PHY.VARC_R,
    type="ELGA",
    location="RIGI",
    components=("TEMP", "HYDR", "SECH", "IRRA", "CORR", "PTOT", "NEUT[2]"),
)

EVARCNR = LocatedComponents(
    phys=PHY.VARC_R,
    type="ELNO",
    components=("TEMP", "HYDR", "SECH", "IRRA", "CORR", "PTOT", "NEUT[2]"),
)


G27NEUTR = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[27]",))

# store all LocatedComponents objects
MODES = objects_from_context(globals(), LocatedComponents)
