from enum import Enum


class TranslationStatus(str, Enum):
    CANDIDATE = "CANDIDATE"
    OFFICIAL = "OFFICIAL"

    def __str__(self) -> str:
        return str(self.value)
