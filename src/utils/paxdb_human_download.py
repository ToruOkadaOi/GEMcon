from bs4 import BeautifulSoup
from urllib.parse import urljoin
import requests
from requests.exceptions import HTTPError
import time
import aiohttp
import asyncio
from aiohttp import ClientTimeout
import pandas as pd
import os
import logging
from rich import print
from rich.console import Console
from rich.table import Table
import yaml
from rich.progress import (
    Progress, BarColumn, TimeElapsedColumn,
    TimeRemainingColumn, DownloadColumn, TransferSpeedColumn
)
timeout = ClientTimeout(total=None)
headers = {"User-Agent": "Mozilla/5.0"}

logger = logging.getLogger("paxdb")
logging.basicConfig(level=logging.INFO)

def get_links():
    base = "https://pax-db.org"
    listing_url = f"{base}/downloads/3.0/datasets/9606/"

    r = requests.get(listing_url)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "html.parser")

    links = []
    for a in soup.find_all("a", href=True):
        href = a["href"]
        if href.endswith(".txt"):
            links.append(urljoin(listing_url, href))

    return links

# download func; In 1 mb chunks
async def download_file(url, filename, chunk_size=1024*1024, retries=3):
    for i in range(retries):
        try:
            connector = aiohttp.TCPConnector(limit=10)
            async with aiohttp.ClientSession(timeout=timeout, headers=headers, connector=connector) as session:
                async with session.get(url) as response:

                    if response.status != 200:
                        logger.error("HTTP error", extra={"status": response.status, "url": url})
                        raise Exception(f"Failed: {response.status}")
                    
                    total = response.content_length

                    from rich.progress import (
                        Progress, BarColumn, TimeElapsedColumn,
                        TimeRemainingColumn, DownloadColumn, TransferSpeedColumn
                    )

                    with Progress(
                        "[progress.description]{task.description}",
                        DownloadColumn(),
                        BarColumn(),
                        TransferSpeedColumn(),
                        TimeRemainingColumn(),
                        TimeElapsedColumn(),
                    ) as progress:

                        task = progress.add_task(
                            f"Downloading {os.path.basename(filename)}",
                            total=total
                        )

                        with open(filename, "wb") as f:
                            async for chunk in response.content.iter_chunked(chunk_size):
                                try:
                                    f.write(chunk)
                                    #raise RuntimeError("test")
                                    progress.update(task, advance=len(chunk))
                                except Exception:
                                    logger.critical("failed to write chunks while downloading", exc_info=True)
                                    raise

            print("Download complete:", filename)
            return

        except Exception as e:
            print(f"Retry {i+1} failed:", e)
            await asyncio.sleep(2 ** i)
save_dir = "data/data_raw/paxdb_human_data" # Download from proj. root
os.makedirs(save_dir, exist_ok=True)

async def download_all(urls):
    tasks = [] # to store the coroutines
    for url in urls:
        filename = os.path.join(save_dir, os.path.basename(url))
        tasks.append(download_file(url, filename)) # use the previous download func.
    
    await asyncio.gather(*tasks) # unpack cor. list

if __name__ == "__main__":
    urls = get_links()
    asyncio.run(download_all(urls))