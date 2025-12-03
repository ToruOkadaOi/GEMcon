from http import HTTPStatus
from io import BytesIO
from typing import Any

import httpx

from ... import errors
from ...client import AuthenticatedClient, Client
from ...types import File, Response


def _get_kwargs(
    id: str,
    type_: str,
) -> dict[str, Any]:
    _kwargs: dict[str, Any] = {
        "method": "get",
        "url": f"/annotation/{id}/download/{type_}",
    }

    return _kwargs


def _parse_response(*, client: AuthenticatedClient | Client, response: httpx.Response) -> File | None:
    if response.status_code == 200:
        response_200 = File(payload=BytesIO(response.json()))

        return response_200

    if client.raise_on_unexpected_status:
        raise errors.UnexpectedStatus(response.status_code, response.content)
    else:
        return None


def _build_response(*, client: AuthenticatedClient | Client, response: httpx.Response) -> Response[File]:
    return Response(
        status_code=HTTPStatus(response.status_code),
        content=response.content,
        headers=response.headers,
        parsed=_parse_response(client=client, response=response),
    )


def sync_detailed(
    id: str,
    type_: str,
    *,
    client: AuthenticatedClient | Client,
) -> Response[File]:
    """This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

     This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

    Args:
        id (str):
        type_ (str):

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        Response[File]
    """

    kwargs = _get_kwargs(
        id=id,
        type_=type_,
    )

    response = client.get_httpx_client().request(
        **kwargs,
    )

    return _build_response(client=client, response=response)


def sync(
    id: str,
    type_: str,
    *,
    client: AuthenticatedClient | Client,
) -> File | None:
    """This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

     This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

    Args:
        id (str):
        type_ (str):

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        File
    """

    return sync_detailed(
        id=id,
        type_=type_,
        client=client,
    ).parsed


async def asyncio_detailed(
    id: str,
    type_: str,
    *,
    client: AuthenticatedClient | Client,
) -> Response[File]:
    """This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

     This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

    Args:
        id (str):
        type_ (str):

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        Response[File]
    """

    kwargs = _get_kwargs(
        id=id,
        type_=type_,
    )

    response = await client.get_async_httpx_client().request(**kwargs)

    return _build_response(client=client, response=response)


async def asyncio(
    id: str,
    type_: str,
    *,
    client: AuthenticatedClient | Client,
) -> File | None:
    """This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

     This is our base controller for annotations that deals with different ontology term types and
    returns a defined annotation schema.

    Args:
        id (str):
        type_ (str):

    Raises:
        errors.UnexpectedStatus: If the server returns an undocumented status code and Client.raise_on_unexpected_status is True.
        httpx.TimeoutException: If the request takes longer than Client.timeout.

    Returns:
        File
    """

    return (
        await asyncio_detailed(
            id=id,
            type_=type_,
            client=client,
        )
    ).parsed
