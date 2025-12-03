from http import HTTPStatus
from typing import Any

import httpx

from ... import errors
from ...client import AuthenticatedClient, Client
from ...models.search_entity_response_200 import SearchEntityResponse200
from ...types import UNSET, Response


def _get_kwargs(
    entity: str,
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
        "url": f"/search/{entity}",
        "params": params,
    }

    return _kwargs


def _parse_response(
    *, client: AuthenticatedClient | Client, response: httpx.Response
) -> SearchEntityResponse200 | None:
    if response.status_code == 200:
        response_200 = SearchEntityResponse200.from_dict(response.json())

        return response_200

    if client.raise_on_unexpected_status:
        raise errors.UnexpectedStatus(response.status_code, response.content)
    else:
        return None


def _build_response(
    *, client: AuthenticatedClient | Client, response: httpx.Response
) -> Response[SearchEntityResponse200]:
    return Response(
        status_code=HTTPStatus(response.status_code),
        content=response.content,
        headers=response.headers,
        parsed=_parse_response(client=client, response=response),
    )


def sync_detailed(
    entity: str,
    *,
    client: AuthenticatedClient | Client,
    q: str,
    page: float = 0.0,
    limit: float = 10.0,
) -> Response[SearchEntityResponse200]:
    """
    Args:
        entity (str):
        q (str):
        page (float):  Default: 0.0.
        limit (float):  Default: 10.0.

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        Response[SearchEntityResponse200]
    """

    kwargs = _get_kwargs(
        entity=entity,
        q=q,
        page=page,
        limit=limit,
    )

    response = client.get_httpx_client().request(
        **kwargs,
    )

    return _build_response(client=client, response=response)


def sync(
    entity: str,
    *,
    client: AuthenticatedClient | Client,
    q: str,
    page: float = 0.0,
    limit: float = 10.0,
) -> SearchEntityResponse200 | None:
    """
    Args:
        entity (str):
        q (str):
        page (float):  Default: 0.0.
        limit (float):  Default: 10.0.

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        SearchEntityResponse200
    """

    return sync_detailed(
        entity=entity,
        client=client,
        q=q,
        page=page,
        limit=limit,
    ).parsed


async def asyncio_detailed(
    entity: str,
    *,
    client: AuthenticatedClient | Client,
    q: str,
    page: float = 0.0,
    limit: float = 10.0,
) -> Response[SearchEntityResponse200]:
    """
    Args:
        entity (str):
        q (str):
        page (float):  Default: 0.0.
        limit (float):  Default: 10.0.

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        Response[SearchEntityResponse200]
    """

    kwargs = _get_kwargs(
        entity=entity,
        q=q,
        page=page,
        limit=limit,
    )

    response = await client.get_async_httpx_client().request(**kwargs)

    return _build_response(client=client, response=response)


async def asyncio(
    entity: str,
    *,
    client: AuthenticatedClient | Client,
    q: str,
    page: float = 0.0,
    limit: float = 10.0,
) -> SearchEntityResponse200 | None:
    """
    Args:
        entity (str):
        q (str):
        page (float):  Default: 0.0.
        limit (float):  Default: 10.0.

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        SearchEntityResponse200
    """

    return (
        await asyncio_detailed(
            entity=entity,
            client=client,
            q=q,
            page=page,
            limit=limit,
        )
    ).parsed
