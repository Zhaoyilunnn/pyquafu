from iqcs.utils.platforms import get_homedir
import os
from .exceptions import UserError
import requests
from .client  import *
import time
import hashlib

def save_token(apitoken):
    """
    Save your apitoken associate your Iqcs account.
    """
    homedir = get_homedir()
    file_dir = homedir + "/.iqcs/"
    if not os.path.exists(file_dir):
        os.mkdir(file_dir)
    with open(file_dir + "api", "w") as f:
        f.write(apitoken+"\n")

def gen_token():
    key = 'ZGFiMjlk=ODY6NTliMz3YWRhNTU0MzJiMT'
    thisH = time.strftime('%Y%m%d%H',time.localtime(time.time()))
    m = hashlib.md5()   
    m.update((thisH+key).encode('utf-8'))   
    md5Str = m.hexdigest()
    r = requests.post(gentoken_url , data = {"key":md5Str})
    res = json.loads(r.content)
    return res['token']
    
def load_token():
    """
    Load Iqcs account.
    """
    homedir = get_homedir()
    file_dir = homedir + "/.iqcs/"
    try:
        with open(file_dir + "api", "r") as f:
            data = f.readlines()
            token = data[0].strip("\n")
            return token
    except:
        raise UserError("User configure error. Please set up your token.")
    