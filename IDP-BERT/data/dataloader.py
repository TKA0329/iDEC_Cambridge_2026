from torch.utils.data import DataLoader
from transformers import BertTokenizerFast
import numpy as np
from data.dataset import IDPBERTDataset


def select_label(train_y, val_y, test_y, label):
    if label == 'rog':
        train_y = train_y[:, 0]
        val_y = val_y[:, 0]
        test_y = test_y[:, 0]
    elif label == 'cv':
        train_y = train_y[:, 1]
        val_y = val_y[:, 1]
        test_y = test_y[:, 1]
    elif label == 'tau':
        train_y = train_y[:, 2]
        val_y = val_y[:, 2]
        test_y = test_y[:, 2]

    return train_y, val_y, test_y


def load_data(config):
    with np.load(f'./data/{config["split"]}/train.npz') as train,\
         np.load(f'./data/{config["split"]}/val.npz') as val,\
         np.load(f'./data/{config["split"]}/test.npz') as test:
        train_X, train_y = train['X'], train['y']
        val_X, val_y = val['X'], val['y']
        test_X, test_y = test['X'], test['y']

    train_y, val_y, test_y = select_label(train_y, val_y, test_y, config['label'])

    tokenizer = BertTokenizerFast.from_pretrained('Rostlab/prot_bert_bfd', do_lower_case=False)

    encoded_tr = tokenizer(list(train_X), padding='max_length', max_length=602, return_attention_mask=True)
    train_encodings, train_attention_masks = encoded_tr['input_ids'], encoded_tr['attention_mask']

    encoded_vl = tokenizer(list(val_X), padding='max_length', max_length=602, return_attention_mask=True)
    val_encodings, val_attention_masks = encoded_vl['input_ids'], encoded_vl['attention_mask']

    encoded_te = tokenizer(list(test_X), padding='max_length', max_length=602, return_attention_mask=True)
    test_encodings, test_attention_masks = encoded_te['input_ids'], encoded_te['attention_mask']

    train_dataset = IDPBERTDataset(input_ids=train_encodings, attention_masks=train_attention_masks, labels=train_y)
    val_dataset = IDPBERTDataset(input_ids=val_encodings, attention_masks=val_attention_masks, labels=val_y)
    test_dataset = IDPBERTDataset(input_ids=test_encodings, attention_masks=test_attention_masks, labels=test_y)

    train_dataloader = DataLoader(
        train_dataset,
        batch_size=config['batch_size'],
        shuffle=True
    )

    val_dataloader = DataLoader(
        val_dataset,
        batch_size=config['batch_size'],
        shuffle=False
    )

    test_dataloader = DataLoader(
        test_dataset,
        batch_size=config['batch_size'],
        shuffle=False
    )

    print('Batch size: ', config['batch_size'])

    print('Train dataset samples: ', len(train_dataset))
    print('Validation dataset samples: ', len(val_dataset))
    print('Test dataset samples: ', len(test_dataset))

    print('Train dataset batches: ', len(train_dataloader))
    print('Validation dataset batches: ', len(val_dataloader))
    print('Test dataset batches: ', len(test_dataloader))

    print()

    return train_dataloader, val_dataloader, test_dataloader
