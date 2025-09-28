import os
import re
import httpx
import typer


class MetadataService:
    @staticmethod
    def download_metadata_jsons(BASE_URL: str, OUT_DIR: str):
        os.makedirs(OUT_DIR, exist_ok=True)

        with httpx.Client(follow_redirects=True, timeout=60.0) as client:
            resp = client.get(BASE_URL)
            resp.raise_for_status()
            html = resp.text

            json_links = re.findall(r'href="([^"]+\.json)"', html)
            typer.echo(f"Found {len(json_links)} JSON files")

            for link in json_links:
                url = BASE_URL + link
                local_path = os.path.join(OUT_DIR, link)

                if os.path.exists(local_path):
                    typer.echo(f"Skipping {link}, already exists")
                    continue

                typer.echo(f"Downloading {url}")
                r = client.get(url)
                r.raise_for_status()
                with open(local_path, "wb") as f:
                    f.write(r.content)

        typer.echo("✅ All JSON files downloaded")

if __name__ == "__main__":
    BASE_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/"
    OUT_DIR = "metadata_jsons"

    MetadataService.download_metadata_jsons(BASE_URL, OUT_DIR)