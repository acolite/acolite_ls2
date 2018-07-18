## def download_file
## download_file with authorisation option
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-14
## modifications:

def download_file(url, file, auth=None, session=None):

    import threading
    import requests
    import time
 
    import os
    file_path = os.path.abspath(file)

    start = time.time()

    def get_url(url, auth=None, stream=True):
        r = requests.get(url, stream=True, auth=auth)

        if (r.ok):
            with open(file_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=1024): 
                    if chunk: # filter out keep-alive new chunks
                        f.write(chunk)
        else:
            print(r.text)
            raise Exception("File download failed")

    threads = [threading.Thread(target=get_url, args=(url,))]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    print("Elapsed Time: {}".format(time.time() - start))
