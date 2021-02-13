FROM continuumio/anaconda3

WORKDIR /bayespairing

COPY . .

RUN chmod +x boot.sh
RUN conda env create -f environment.yml
SHELL ["conda", "run", "-n", "bayespairing", "/bin/bash", "-c"]

EXPOSE 5000

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bayespairing", "./boot.sh"]