import torch
from tqdm import tqdm
from sklearn.metrics import r2_score
from transformers import BertTokenizerFast


def train(model, dataloader, optimizer, criterion, device):
    model.train()
    total_loss = 0.0

    ground_truth = []
    predictions = []

    for batch in tqdm(dataloader):
        inputs = batch['input_ids'].to(device)
        attention_mask = batch['attention_mask'].to(device)
        labels = batch['labels'].to(device)

        optimizer.zero_grad()

        outputs = model(inputs, attention_mask)
        loss = criterion(outputs.squeeze(1), labels)

        loss.backward()
        optimizer.step()

        total_loss += loss.item()

        predictions.extend(outputs.detach().cpu().tolist())
        ground_truth.extend(labels.detach().cpu().tolist())

    return total_loss / len(dataloader), r2_score(ground_truth, predictions)


def validate(model, dataloader, criterion, device):
    model.eval()
    total_loss = 0.0

    ground_truth = []
    predictions = []

    for batch in tqdm(dataloader):
        inputs = batch['input_ids'].to(device)
        attention_mask = batch['attention_mask'].to(device)
        labels = batch['labels'].to(device)

        with torch.inference_mode():
            outputs = model(inputs, attention_mask)
            loss = criterion(outputs.squeeze(1), labels)

        total_loss += loss.item()

        predictions.extend(outputs.detach().cpu().tolist())
        ground_truth.extend(labels.detach().cpu().tolist())

    return total_loss / len(dataloader), r2_score(ground_truth, predictions)


def test(model, dataloader, device):
    model.eval()

    ground_truth = []
    predictions = []

    for batch in tqdm(dataloader):
        inputs = batch['input_ids'].to(device)
        attention_mask = batch['attention_mask'].to(device)
        labels = batch['labels']

        with torch.inference_mode():
            outputs = model(inputs, attention_mask).squeeze()

        if outputs.dim()==0:
            outputs = outputs.unsqueeze(0)
        if labels.dim()==0:
            labels = labels.unsqueeze(0)

        predictions.extend(outputs.detach().cpu().tolist())
        ground_truth.extend(labels.detach().cpu().tolist())

    return r2_score(ground_truth, predictions)

def get_predictions(model, dataloader, device):
    model.eval()
    predictions = []

    for batch in tqdm(dataloader):
        inputs = batch['input_ids'].to(device)
        attention_mask = batch['attention_mask'].to(device)

        with torch.inference_mode():
            outputs = model(inputs, attention_mask).squeeze()

        if outputs.dim()==0:
            outputs = outputs.unsqueeze(0)

        predictions.extend(outputs.detach().cpu().tolist())

    return predictions
