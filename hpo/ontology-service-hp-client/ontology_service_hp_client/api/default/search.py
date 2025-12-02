from http import HTTPStatus
from typing import Any

import httpx

from ... import errors
from ...client import AuthenticatedClient, Client
from ...models.search_dto import SearchDto
from ...types import UNSET, Response


def _get_kwargs(
    *,
    q: str,
    page: float = 0.0,
    limit: float = 10.0,
) -> dict[str, Any]:
    params: dict[str, Any] = {}

    params["q"] = q

    params["page"] = page

    params["limit"] = limit

    params = {k: v for k, v in params.items() if v is not UNSET and v is not None}

    _kwargs: dict[str, Any] = {
        "method": "get",
        "url": "/api/hp/search",
        "params": params,
    }

    return _kwargs


def _parse_response(*, client: AuthenticatedClient | Client, response: httpx.Response) -> SearchDto | None:
    if response.status_code == 200:
        response_200 = SearchDto.from_dict(response.json())

        return response_200

    if client.raise_on_unexpected_status:
        raise errors.UnexpectedStatus(response.status_code, response.content)
    else:
        return None


def _build_response(*, client: AuthenticatedClient | Client, response: httpx.Response) -> Response[SearchDto]:
    return Response(
        status_code=HTTPStatus(response.status_code),
        content=response.content,
        headers=response.headers,
        parsed=_parse_response(client=client, response=response),
    )


def sync_detailed(
    *,
    client: AuthenticatedClient | Client,
    q: str,
    page: float = 0.0,
    limit: float = 10.0,
) -> Response[SearchDto]:
    """Search for ontology term, id, synonym

     Search for ontology term, id, synonym

    Args:
        q (str):
        page (float):  Default: 0.0.
        limit (float):  Default: 10.0.

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        Response[SearchDto]
    """

    kwargs = _get_kwargs(
        q=q,
        page=page,
        limit=limit,
    )

    response = client.get_httpx_client().request(
        **kwargs,
    )

    return _build_response(client=client, response=response)


def sync(
    *,
    client: AuthenticatedClient | Client,
    q: str,
    page: float = 0.0,
    limit: float = 10.0,
) -> SearchDto | None:
    """Search for ontology term, id, synonym

     Search for ontology term, id, synonym

    Args:
        q (str):
        page (float):  Default: 0.0.
        limit (float):  Default: 10.0.

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        SearchDto
    """

    return sync_detailed(
        client=client,
        q=q,
        page=page,
        limit=limit,
    ).parsed


async def asyncio_detailed(
    *,
    client: AuthenticatedClient | Client,
    q: str,
    page: float = 0.0,
    limit: float = 10.0,
) -> Response[SearchDto]:
    """Search for ontology term, id, synonym

     Search for ontology term, id, synonym

    Args:
        q (str):
        page (float):  Default: 0.0.
        limit (float):  Default: 10.0.

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        Response[SearchDto]
    """

    kwargs = _get_kwargs(
        q=q,
        page=page,
        limit=limit,
    )

    response = await client.get_async_httpx_client().request(**kwargs)

    return _build_response(client=client, response=response)


async def asyncio(
    *,
    client: AuthenticatedClient | Client,
    q: str,
    page: float = 0.0,
    limit: float = 10.0,
) -> SearchDto | None:
    """Search for ontology term, id, synonym

     Search for ontology term, id, synonym

    Args:
        q (str):
        page (float):  Default: 0.0.
        limit (float):  Default: 10.0.

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        SearchDto
    """

    return (
        await asyncio_detailed(
            client=client,
            q=q,
            page=page,
            limit=limit,
        )
    ).parsed
