import logging
import math
import os
from collections import Counter
from dataclasses import dataclass
from enum import Enum
from typing import List, Optional

import numpy.typing

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


class Ion(Enum):
    Ag = "Ag"
    Au = "Au"
    Ba = "Ba"
    Ca = "Ca"
    Co = "Co"
    Cs = "Cs"
    Cu = "Cu"
    Eu = "Eu"
    Fe = "Fe"
    Ir = "Ir"
    K = "K"
    Li = "Li"
    Mg = "Mg"
    Mn = "Mn"
    Na = "Na"
    Ni = "Ni"
    Os = "Os"
    Pb = "Pb"
    Pt = "Pt"
    Ru = "Ru"
    Sr = "Sr"
    Tl = "Tl"
    V = "V"
    Zn = "Zn"


class ONZ(Enum):
    O_PLUS = "O+"
    O_MINUS = "O-"
    N_PLUS = "N+"
    N_MINUS = "N-"
    Z_PLUS = "Z+"
    Z_MINUS = "Z-"

    def score(self):
        if self == ONZ.O_PLUS:
            return 0
        elif self == ONZ.O_MINUS:
            return 1
        elif self == ONZ.N_PLUS:
            return 2
        elif self == ONZ.N_MINUS:
            return 3
        elif self == ONZ.Z_PLUS:
            return 4
        elif self == ONZ.Z_MINUS:
            return 5
        return 6


class ONZM(Enum):
    OP_PLUS = "Op+"
    OP_MINUS = "Op-"
    OP_STAR = "Op*"
    OA_PLUS = "Oa+"
    OA_MINUS = "Oa-"
    OA_STAR = "Oa*"
    OH_PLUS = "Oh+"
    OH_MINUS = "Oh-"
    OH_STAR = "Oh*"
    NP_PLUS = "Np+"
    NP_MINUS = "Np-"
    NP_STAR = "Np*"
    NA_PLUS = "Na+"
    NA_MINUS = "Na-"
    NA_STAR = "Na*"
    NH_PLUS = "Nh+"
    NH_MINUS = "Nh-"
    NH_STAR = "Nh*"
    ZP_PLUS = "Zp+"
    ZP_MINUS = "Zp-"
    ZP_STAR = "Zp*"
    ZA_PLUS = "Za+"
    ZA_MINUS = "Za-"
    ZA_STAR = "Za*"
    ZH_PLUS = "Zh+"
    ZH_MINUS = "Zh-"
    ZH_STAR = "Zh*"
    MP_PLUS = "Mp+"
    MP_MINUS = "Mp-"
    MP_STAR = "Mp*"
    MA_PLUS = "Ma+"
    MA_MINUS = "Ma-"
    MA_STAR = "Ma*"
    MH_PLUS = "Mh+"
    MH_MINUS = "Mh-"
    MH_STAR = "Mh*"

    @staticmethod
    def from_value(value: str):
        for candidate in ONZM:
            if candidate.value == value:
                return candidate
        raise RuntimeError(f"Failed to match {value} to an ONZM class")


class GbaTetradClassification(Enum):
    Ia = "Ia"
    IIa = "IIa"
    IIIa = "IIIa"
    IVa = "IVa"
    Va = "Va"
    VIa = "VIa"
    VIIa = "VIIa"
    VIIIa = "VIIIa"
    Ib = "Ib"
    IIb = "IIb"
    IIIb = "IIIb"
    IVb = "IVb"
    Vb = "Vb"
    VIb = "VIb"
    VIIb = "VIIb"
    VIIIb = "VIIIb"


class GbaQuadruplexClassification(Enum):
    I = "I"
    II = "II"
    III = "III"
    IV = "IV"
    V = "V"
    VI = "VI"
    VII = "VII"
    VIII = "VIII"


class LoopClassification(Enum):
    _1a = "1a"
    _2a = "2a"
    _3a = "3a"
    _4a = "4a"
    _5a = "5a"
    _6a = "6a"
    _7a = "7a"
    _8a = "8a"
    _9a = "9a"
    _10a = "10a"
    _11a = "11a"
    _12a = "12a"
    _13a = "13a"
    _1b = "1b"
    _2b = "2b"
    _3b = "3b"
    _4b = "4b"
    _5b = "5b"
    _6b = "6b"
    _7b = "7b"
    _8b = "8b"
    _9b = "9b"
    _10b = "10b"
    _11b = "11b"
    _12b = "12b"
    _13b = "13b"
    invalid = "n/a"

    def loop_progression(self) -> str:
        if self == LoopClassification._1a:
            return "-(ppp)"
        elif self == LoopClassification._1b:
            return "+(ppp)"
        elif self == LoopClassification._2a:
            return "-(ppl)"
        elif self == LoopClassification._2b:
            return "+(ppl)"
        elif self == LoopClassification._3a:
            return "-(plp)"
        elif self == LoopClassification._3b:
            return "+(plp)"
        elif self == LoopClassification._4a:
            return "-(lpp)"
        elif self == LoopClassification._4b:
            return "+(lpp)"
        elif self == LoopClassification._5a:
            return "(-pd+p)"
        elif self == LoopClassification._5b:
            return "(+pd-p)"
        elif self == LoopClassification._6a:
            return "-(lll)"
        elif self == LoopClassification._6b:
            return "+(lll)"
        elif self == LoopClassification._7a:
            return "-(llp)"
        elif self == LoopClassification._7b:
            return "+(llp)"
        elif self == LoopClassification._8a:
            return "-(lpl)"
        elif self == LoopClassification._8b:
            return "+(lpl)"
        elif self == LoopClassification._9a:
            return "-(pll)"
        elif self == LoopClassification._9b:
            return "+(pll)"
        elif self == LoopClassification._10a:
            return "(-pd+l)"
        elif self == LoopClassification._10b:
            return "(+pd-l)"
        elif self == LoopClassification._11a:
            return "(-ld+l)"
        elif self == LoopClassification._11b:
            return "(+ld-l)"
        elif self == LoopClassification._12a:
            return "(d-pd)"
        elif self == LoopClassification._12b:
            return "(d+pd)"
        elif self == LoopClassification._13a:
            return "(-ld+p)"
        elif self == LoopClassification._13b:
            return "(+ld-p)"
        elif self == LoopClassification.invalid:
            return "n/a"
        raise RuntimeError(f"Failed to get string representation of {self}")

    @staticmethod
    def from_value(value: str):
        for candidate in LoopClassification:
            if candidate.value == value:
                return candidate
        raise RuntimeError(f"Failed to match {value} to a LoopClassification class")


class LoopType(Enum):
    diagonal = "diagonal"
    propeller_plus = "propeller+"
    propeller_minus = "propeller-"
    lateral_plus = "lateral+"
    lateral_minus = "lateral-"

    @staticmethod
    def from_value(value: str):
        for loop_type in LoopType:
            if loop_type.value == value:
                return loop_type
        raise RuntimeError(f"Failed to match {value} to a LoopType class")


class Direction(Enum):
    parallel = "parallel"
    antiparallel = "antiparallel"
    hybrid = "hybrid"
