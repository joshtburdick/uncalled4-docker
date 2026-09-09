# uncalled4-docker
Docker image to run uncalled4

## Building the image

Build the image using the following command:

```bash
docker build -t uncalled4 .
```

## Running the image

To run `uncalled4`, you should mount your local data directory to the `/data` volume inside the container:

```bash
docker run --rm -v /path/to/your/data:/data uncalled4 [command]
```

For example, to see the help menu:

```bash
docker run --rm uncalled4 --help
```
