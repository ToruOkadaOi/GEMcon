from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, TypeVar

from attrs import define as _attrs_define
from attrs import field as _attrs_field

from ..types import UNSET, Unset

if TYPE_CHECKING:
    from ..models.translation import Translation


T = TypeVar("T", bound="SimpleOntologyTerm")


@_attrs_define
class SimpleOntologyTerm:
    """
    Attributes:
        id (str | Unset):
        name (str | Unset):
        translations (list[Translation] | Unset):
        descendant_count (int | Unset):
    """

    id: str | Unset = UNSET
    name: str | Unset = UNSET
    translations: list[Translation] | Unset = UNSET
    descendant_count: int | Unset = UNSET
    additional_properties: dict[str, Any] = _attrs_field(init=False, factory=dict)

    def to_dict(self) -> dict[str, Any]:
        id = self.id

        name = self.name

        translations: list[dict[str, Any]] | Unset = UNSET
        if not isinstance(self.translations, Unset):
            translations = []
            for translations_item_data in self.translations:
                translations_item = translations_item_data.to_dict()
                translations.append(translations_item)

        descendant_count = self.descendant_count

        field_dict: dict[str, Any] = {}
        field_dict.update(self.additional_properties)
        field_dict.update({})
        if id is not UNSET:
            field_dict["id"] = id
        if name is not UNSET:
            field_dict["name"] = name
        if translations is not UNSET:
            field_dict["translations"] = translations
        if descendant_count is not UNSET:
            field_dict["descendantCount"] = descendant_count

        return field_dict

    @classmethod
    def from_dict(cls: type[T], src_dict: Mapping[str, Any]) -> T:
        from ..models.translation import Translation

        d = dict(src_dict)
        id = d.pop("id", UNSET)

        name = d.pop("name", UNSET)

        _translations = d.pop("translations", UNSET)
        translations: list[Translation] | Unset = UNSET
        if _translations is not UNSET:
            translations = []
            for translations_item_data in _translations:
                translations_item = Translation.from_dict(translations_item_data)

                translations.append(translations_item)

        descendant_count = d.pop("descendantCount", UNSET)

        simple_ontology_term = cls(
            id=id,
            name=name,
            translations=translations,
            descendant_count=descendant_count,
        )

        simple_ontology_term.additional_properties = d
        return simple_ontology_term

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
