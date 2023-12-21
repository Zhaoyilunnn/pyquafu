from iqcs.exceptions import CircuitError, ServerError, CompileError
import requests
import json
from iqcs.utils.convertors import circuit_to_json
from typing import List

# base_url = 'http://127.0.0.1:9999/'
base_url = 'http://10.58.96.25:9999/'
submit_url = base_url + 'experiments/task/submit'
query_url = base_url + 'experiments/status/query'
info_url = base_url + 'machines/info/query'
gentoken_url = base_url + 'token/generate'
validate_url = base_url + 'experiment/validate'


def get_bakcend_info():
     r = requests.post(info_url)
     print(str(r.content, encoding = "utf-8"))