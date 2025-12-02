from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, TypeVar

from attrs import define as _attrs_define
from attrs import field as _attrs_field

if TYPE_CHECKING:
    from ..models.ontology_term import OntologyTerm


T = TypeVar("T", bound="SearchDto")


@_attrs_define
class SearchDto:
    """
    Attributes:
        terms (list[OntologyTerm]):
        total_count (int):
    """

    terms: list[OntologyTerm]
    total_count: int
    additional_properties: dict[str, Any] = _attrs_field(init=False, factory=dict)

    def to_dict(self) -> dict[str, Any]:
        terms = []
        for terms_item_data in self.terms:
            terms_item = terms_item_data.to_dict()
            terms.append(terms_item)

        total_count = self.total_count

        field_dict: dict[str, Any] = {}
        field_dict.update(self.additional_properties)
        field_dict.update(
            {
                "terms": terms,
                "totalCount": total_count,
            }
        )

        return field_dict

    @classmethod
    def from_dict(cls: type[T], src_dict: Mapping[str, Any]) -> T:
        from ..models.ontology_term import OntologyTerm

        d = dict(src_dict)
        terms = []
        _terms = d.pop("terms")
        for terms_item_data in _terms:
            terms_item = OntologyTerm.from_dict(terms_item_data)

            terms.append(terms_item)

        total_count = d.pop("totalCount")

        search_dto = cls(
            terms=terms,
            total_count=total_count,
        )

        search_dto.additional_properties = d
        return search_dto

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
