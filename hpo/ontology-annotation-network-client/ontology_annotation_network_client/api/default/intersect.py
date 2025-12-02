from http import HTTPStatus
from typing import Any

import httpx

from ... import errors
from ...client import AuthenticatedClient, Client
from ...models.intersect_response_200 import IntersectResponse200
from ...types import UNSET, Response


def _get_kwargs(
    entity: str,
    *,
    p: str,
) -> dict[str, Any]:
    params: dict[str, Any] = {}

    params["p"] = p

    params = {k: v for k, v in params.items() if v is not UNSET and v is not None}

    _kwargs: dict[str, Any] = {
        "method": "get",
        "url": f"/search/{entity}/intersect",
        "params": params,
    }

    return _kwargs


def _parse_response(*, client: AuthenticatedClient | Client, response: httpx.Response) -> IntersectResponse200 | None:
    if response.status_code == 200:
        response_200 = IntersectResponse200.from_dict(response.json())

        return response_200

    if client.raise_on_unexpected_status:
        raise errors.UnexpectedStatus(response.status_code, response.content)
    else:
        return None


def _build_response(
    *, client: AuthenticatedClient | Client, response: httpx.Response
) -> Response[IntersectResponse200]:
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
    p: str,
) -> Response[IntersectResponse200]:
    """This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

     This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

    Args:
        entity (str):
        p (str):

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        Response[IntersectResponse200]
    """

    kwargs = _get_kwargs(
        entity=entity,
        p=p,
    )

    response = client.get_httpx_client().request(
        **kwargs,
    )

    return _build_response(client=client, response=response)


def sync(
    entity: str,
    *,
    client: AuthenticatedClient | Client,
    p: str,
) -> IntersectResponse200 | None:
    """This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

     This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

    Args:
        entity (str):
        p (str):

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        IntersectResponse200
    """

    return sync_detailed(
        entity=entity,
        client=client,
        p=p,
    ).parsed


async def asyncio_detailed(
    entity: str,
    *,
    client: AuthenticatedClient | Client,
    p: str,
) -> Response[IntersectResponse200]:
    """This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

     This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

    Args:
        entity (str):
        p (str):

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        Response[IntersectResponse200]
    """

    kwargs = _get_kwargs(
        entity=entity,
        p=p,
    )

    response = await client.get_async_httpx_client().request(**kwargs)

    return _build_response(client=client, response=response)


async def asyncio(
    entity: str,
    *,
    client: AuthenticatedClient | Client,
    p: str,
) -> IntersectResponse200 | None:
    """This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

     This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

    Args:
        entity (str):
        p (str):

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        IntersectResponse200
    """

    return (
        await asyncio_detailed(
            entity=entity,
            client=client,
            p=p,
        )
    ).parsed
