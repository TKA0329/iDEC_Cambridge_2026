import torch
import yaml
import importlib.util
import os
from transformers import BertTokenizer, BertModel
from transformers import BertConfig

import random
from huggingface_hub import hf_hub_download
import numpy as np
import csv
import re
import pickle
from sklearn.preprocessing import PowerTransformer

import pandas as pd 
from tqdm import tqdm
import numpy as np

# numpy compatibility + PyTorch safe loading 
# the location changed between NumPy versions
try:
    from numpy._core.multiarray import scalar  # numpy 2.x
except ImportError:
    from numpy.core.multiarray import scalar   # numpy 1.x

# PyTorch's safer model deserialisation 
# if .pt file contains something involving a NumPy scalar, 
# this tells PyTorch this particular class/function is allowed when unpickling this model
torch.serialization.add_safe_globals([scalar])


def parse_csv(filename) -> list:
    aa_variant = []
    with open(f"{filename}", "r", encoding = "utf-8") as file: #open the csv file
        reader = csv.DictReader(file)
        for row in reader:
            aa = row["aa_variant"]
            aa_variant.append(aa)
    return aa_variant

with open("data/data.pkl", "rb") as f:   # only load if you trust the source
    raw = pickle.load(f)

scalers = {}
for k in ["rog", "cv", "tau"]:
    s = PowerTransformer()
    s.fit(raw[k].reshape(-1, 1))
    scalers[k] = s

def to_real(name, value):
    return scalers[name].inverse_transform([[value]])[0, 0]


# making experiments reproducible 
# run the same computation as consistently as possible between runs 
def seed_everything(seed=42):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    # telling PyTorch / cuDNN to prefer deterministic algorithms
    torch.backends.cudnn.deterministic = True 
    # turns of cuDNN automatic benchmarking that can choose diff algo depending on the hardware/input
    torch.backends.cudnn.benchmark = False

# what models im using 
MODELS = ["cv", "rog", "tau"]
BASE_DIR = "best_models"

def load_models():
    models = {}
# looks like:
#     # models = {
#     "cv": <CV model>,
#     "rog": <ROG model>,
#     "tau": <TAU model>
# }

# loop over the three models 
    for model_dir in MODELS:
        # where the model weights should be 
        local_path = os.path.join(BASE_DIR, model_dir, "model.pt")
        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        # downloading the model if nucessary 
        if not os.path.exists(local_path):
            hf_hub_download(
                repo_id="dsadasiv/IDP-BERT",
                filename=f"best_models/{model_dir}/model.pt",
                local_dir=".",
                local_dir_use_symlinks=False
            )

        # contains model-specific configuration 
        config_path = os.path.join(model_dir, "config.yaml")
        # contains the trained parameters 
        weights_path = os.path.join(BASE_DIR, model_dir, "model.pt")
        # contains the actual Python definition of the IDP-Bert neural network 
        network_path = os.path.join(model_dir, "network.py")

        # Take this python file and make it into a module 
        spec = importlib.util.spec_from_file_location(f"network_{model_dir}", network_path)
        net_module = importlib.util.module_from_spec(spec) # creates the module object 
        spec.loader.exec_module(net_module) # executes network.py 

        with open(config_path, 'r') as f:
            idp_config = yaml.safe_load(f) # reads params and convert it into a Python dict 

        # The architecture of the bert component 
        # telling HF to construct a BERT with this architecture 
        bert_config = BertConfig(
                vocab_size=30522,
                hidden_size=256,
                num_hidden_layers=16,
                intermediate_size=3072,
                max_position_embeddings=602,
                type_vocab_size=2,
                num_attention_heads=16,   # from repo's config 
            )
        model = net_module.IDPBERT(bert_config, idp_config, get_embeddings=False) # creates the actual model 
        # creates the architecture but does NOT YET contain the trained weights 
        # at this point its an untrained/initialised model with the correct structure 

        ckpt = torch.load(weights_path, map_location="cpu", weights_only=False) # reads model.pt 
        
        

        # mapping parameter names to tensors 
        state_dict = ckpt["model_state_dict"] # takes the trained parameters out of the check point 

        # put those weights into the model / arcitecture 
        # strict = True means parameter names and shapes must match exactly 
        model.load_state_dict(state_dict, strict=True)
        # missing, unexpected = model.load_state_dict(state_dict, strict=False)
        # print(model_dir, "missing:", missing, "unexpected:", unexpected)
        
        # switching to evaluation mode 
        model.eval()

        # storing the model 
        models[model_dir] = model
    
        print(f"{model_dir}: val R^2 = {ckpt['val_r2']}")
    print({k: id(v) for k, v in models.items()})
    
    return models


def device():
    device = "CUDA" if torch.cuda.is_available() else "CPU"
    print(f"**Device:** {device}")
    if torch.cuda.is_available():
        print(f"**GPU:** {torch.cuda.get_device_name(0)}")
    print(f"**PyTorch:** {torch.__version__}")
    print(f"**Models loaded:** {', '.join(MODELS)}")

# loads the tokenizer -> loads/downloads the tokenizer associated with ProtBert 
# converts protein seq from a string into numerical token IDs that Bert understands 
tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert_bfd", do_lower_case=False)


def predict(sequence, model):
    sequence = re.sub(r"[UZOB]", "X", sequence.upper())
    # turns MKT -> M K T
    processed_seq = " ".join(list(sequence.replace(" ", "")))
    # tokenises it 
    inputs = tokenizer(processed_seq, padding='max_length', 
                       max_length=602, 
                       truncation = True,
                       return_attention_mask=True,
                       return_tensors = "pt")
    with torch.no_grad():
        # input_ids -> numerical rep of the aa seq 
        # attention_mask -> tells BERT which positions are actual sequence tokens vs padding 
        prediction = model(inputs['input_ids'], inputs['attention_mask'])
    return prediction.item()

# import numpy as np
# from scipy.stats import pearsonr


# def pkt(models):
#     with open("data/data.pkl", "rb") as f:
#         raw = pickle.load(f)
#         print(type(raw))

#     if isinstance(raw, dict):
#         for k, v in raw.items():
#             shape = getattr(v, "shape", None)
#             length = len(v) if hasattr(v, "__len__") else None
#             print(repr(k), type(v).__name__, "shape:", shape, "len:", length)
#             # peek at the first element
#             try:
#                 print("   first:", repr(v[0])[:100])
#             except Exception:
#                 pass
#     for k, v in raw.items():
#         print(k, type(v), getattr(v, "shape", len(v) if hasattr(v, "__len__") else v))
#         seq_key = "seqs"   # whatever the key actually is
#         n = 50
#         idx = np.random.RandomState(0).choice(len(raw[seq_key]), n, replace=False)

#         for name in ["rog", "cv", "tau"]:
#             pred = np.array([to_real(name, predict(raw[seq_key][i], models[name])) for i in idx])
#             true = np.array([raw[name][i] for i in idx])
#             print(name,
#                 "pearson:", round(pearsonr(pred, true)[0], 3),
#                 "median rel. error:", round(np.median(np.abs(pred - true) / np.abs(true)), 3))

    


def main():
    filename= input("filename: ")
    destination = input("Destination file: ")
    seed_everything(42)
    models = load_models()
    device()
    result_list = []
    aa_variant = parse_csv(filename)
    for sequence in tqdm(aa_variant):
        for name, model in models.items():
            result = predict(sequence, model)
            result_real = to_real(name, result)
            result_list.append({"Model's name": name,
                                "result": result_real,
                                "sequence": sequence})
            # print(f"----------------------------------------")
            # print(f"sequence: {sequence}: ")
            # print(f"Model: {name}; result: {result_real}")
            # print(f"----------------------------------------")
    df = pd.read_csv(filename)
    res = pd.DataFrame(result_list)
    wide = (
    res.pivot_table(index="sequence", columns="Model's name", values="result")
       .rename(columns={"rog": "pred_rog_A", "cv": "pred_cv_JperK", "tau": "pred_tau_fs"}))

    out = df.merge(wide, left_on="aa_variant", right_index=True, how="left")

    # alignment checks
    assert len(out) == len(df)                       # no rows added or lost
    assert out[wide.columns].notna().all().all()     # every row got all 3 predictions

    out.to_csv(f"{destination}", index=False)
    for k in ["rog", "cv", "tau"]:
        print("raw stuff: ")
        print(k, raw[k].min(), raw[k].max())

    # pkt(models)
    print(result_list)
            


if __name__ == "__main__":
    main()