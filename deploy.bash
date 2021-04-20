# stops, removes and builds an image for BP
# runs the container on port 5001

sudo docker container stop bayespairing
sudo docker container rm bayespairing 

sudo docker build -t bayespairing .
sudo docker run -p 5001:5001 -d -t --restart on-failure --network=PipelineNetwork  --name bayespairing bayespairing 
