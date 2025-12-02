from enum import Enum


class Language(str, Enum):
    CS = "CS"
    DE = "DE"
    DTP = "DTP"
    EN = "EN"
    ES = "ES"
    FR = "FR"
    IT = "IT"
    JA = "JA"
    NL = "NL"
    NNA = "NNA"
    PT = "PT"
    TR = "TR"
    TW = "TW"
    ZH = "ZH"

    def __str__(self) -> str:
        return str(self.value)
