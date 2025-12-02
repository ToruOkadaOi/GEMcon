from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, TypeVar, cast

from attrs import define as _attrs_define
from attrs import field as _attrs_field

from ..types import UNSET, Unset

if TYPE_CHECKING:
    from ..models.translation import Translation


T = TypeVar("T", bound="OntologyTerm")


@_attrs_define
class OntologyTerm:
    """
    Attributes:
        descendant_count (int):
        synonyms (list[str]):
        xrefs (list[str]):
        id (str | Unset):
        name (str | Unset):
        definition (str | Unset):
        comment (str | Unset):
        translations (list[Translation] | Unset):
    """

    descendant_count: int
    synonyms: list[str]
    xrefs: list[str]
    id: str | Unset = UNSET
    name: str | Unset = UNSET
    definition: str | Unset = UNSET
    comment: str | Unset = UNSET
    translations: list[Translation] | Unset = UNSET
    additional_properties: dict[str, Any] = _attrs_field(init=False, factory=dict)

    def to_dict(self) -> dict[str, Any]:
        descendant_count = self.descendant_count

        synonyms = self.synonyms

        xrefs = self.xrefs

        id = self.id

        name = self.name

        definition = self.definition

        comment = self.comment

        translations: list[dict[str, Any]] | Unset = UNSET
        if not isinstance(self.translations, Unset):
            translations = []
            for translations_item_data in self.translations:
                translations_item = translations_item_data.to_dict()
                translations.append(translations_item)

        field_dict: dict[str, Any] = {}
        field_dict.update(self.additional_properties)
        field_dict.update(
            {
                "descendantCount": descendant_count,
                "synonyms": synonyms,
                "xrefs": xrefs,
            }
        )
        if id is not UNSET:
            field_dict["id"] = id
        if name is not UNSET:
            field_dict["name"] = name
        if definition is not UNSET:
            field_dict["definition"] = definition
        if comment is not UNSET:
            field_dict["comment"] = comment
        if translations is not UNSET:
            field_dict["translations"] = translations

        return field_dict

    @classmethod
    def from_dict(cls: type[T], src_dict: Mapping[str, Any]) -> T:
        from ..models.translation import Translation

        d = dict(src_dict)
        descendant_count = d.pop("descendantCount")

        synonyms = cast(list[str], d.pop("synonyms"))

        xrefs = cast(list[str], d.pop("xrefs"))

        id = d.pop("id", UNSET)

        name = d.pop("name", UNSET)

        definition = d.pop("definition", UNSET)

        comment = d.pop("comment", UNSET)

        _translations = d.pop("translations", UNSET)
        translations: list[Translation] | Unset = UNSET
        if _translations is not UNSET:
            translations = []
            for translations_item_data in _translations:
                translations_item = Translation.from_dict(translations_item_data)

                translations.append(translations_item)

        ontology_term = cls(
            descendant_count=descendant_count,
            synonyms=synonyms,
            xrefs=xrefs,
            id=id,
            name=name,
            definition=definition,
            comment=comment,
            translations=translations,
        )

        ontology_term.additional_properties = d
        return ontology_term

    @property
    def additional_keys(self) -> list[str]:
        return list(self.additional_properties.keys())

    def __getitem__(self, key: str) -> Any:
        return self.additional_properties[key]

    def __setitem__(self, key: str, value: Any) -> None:
        self.additional_properties[key] = value

    def __delitem__(self, key: str) -> None:
        del self.additional_properties[key]

    def __contains__(self, key: str) -> bool:
        return key in self.additional_properties
