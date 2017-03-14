FROM python:3
ADD ./ /
RUN pip install -r requirements.txt
CMD [ "python", "./upper_quantile_normalize.py","--genefile","genefile","--genmatfile","gen_matrix_file"]