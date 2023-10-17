FROM jupyter/base-notebook:ubuntu-22.04

# Name your environment and choose the python version
ARG env_name=deep-dye-drop
#ARG py_ver=3.8.5
USER root

RUN apt-get update -qq && \
    apt-get install -y git

# You can add additional libraries here
# RUN mamba create --yes -p "${CONDA_DIR}/envs/${env_name}" \
#     python=${py_ver} \
#     'ipykernel' \
#     'jupyterlab' && \
#     mamba clean --all -f -y

# Alternatively, you can comment out the lines above and uncomment those below
# if you'd prefer to use a YAML file present in the docker build context

COPY --chown=${NB_UID}:${NB_GID} env.yaml /tmp/
RUN mamba env create -p "${CONDA_DIR}/envs/${env_name}" -f /tmp/env.yaml && \
    mamba clean --all -f -y



# Create Python kernel and link it to jupyter
RUN "${CONDA_DIR}/envs/${env_name}/bin/python" -m ipykernel install --user --name="${env_name}" && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

# Any additional `pip` installs can be added by using the following line
# Using `mamba` is highly recommended though
# RUN "${CONDA_DIR}/envs/${env_name}/bin/pip" install --no-cache-dir \
#     'flake8'

### install datarail, gr50, and cell_cycle_gating packages from github
RUN "${CONDA_DIR}/envs/${env_name}/bin/pip" install git+https://github.com/datarail/datarail.git@068a7ba35bf36444f708bd5161e4d99ca8266188
RUN "${CONDA_DIR}/envs/${env_name}/bin/pip" install git+https://github.com/datarail/gr_metrics.git@9c181322bb0dbca3459d80840b26104178fdf803#egg=gr50\&subdirectory=SRC/python
RUN "${CONDA_DIR}/envs/${env_name}/bin/pip" install git+https://github.com/datarail/DrugResponse.git@3163393ed6c7b44dc0f6e528c013671bd3dfd249#egg=cell_cycle_gating\&subdirectory=python

# Creating a startup hook, which will activate our custom environment by default in Jupyter Notebook
# More info about startup hooks: https://jupyter-docker-stacks.readthedocs.io/en/latest/using/common.html#startup-hooks
# You can comment this section to keep the default environment in Jupyter Notebook
USER root
RUN activate_custom_env_script=/usr/local/bin/before-notebook.d/activate_custom_env.sh && \
    echo "#!/bin/bash" > ${activate_custom_env_script} && \
    echo "eval \"$(conda shell.bash activate "${env_name}")\"" >> ${activate_custom_env_script} && \
    chmod +x ${activate_custom_env_script}

USER ${NB_UID}

# Making this environment default in Terminal
# You can comment this line to keep the default environment in Terminal
RUN echo "conda activate ${env_name}" >> "${HOME}/.bashrc"
