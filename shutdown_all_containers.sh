#!/bin/bash

# Get a list of running Docker containers
running_containers=$(docker ps -q)

# Check if there are any running containers
if [ -n "$running_containers" ]; then
  echo "Running Docker containers found. Stopping them..."
  
  # Loop through each running container and stop it
  for container_id in $running_containers; do
    docker stop "$container_id"
    docker rm "$container_id"
    echo "Stopped container: $container_id"
  done

else
  echo "No running Docker containers found."
fi