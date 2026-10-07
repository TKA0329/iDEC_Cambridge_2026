import torch
from transformers import BertModel, BertConfig, logging

logging.set_verbosity_error()


class IDPBERT(torch.nn.Module):
    def __init__(self, bert_config, config, get_embeddings):
        super(IDPBERT, self).__init__()

        self.protbert = BertModel.from_pretrained(
            'Rostlab/prot_bert_bfd',
            config=bert_config,
            ignore_mismatched_sizes=True
        )

        layers = []
        for _ in range(config['network']['fc_layers']-1):
            layers.append(torch.nn.Linear(bert_config.hidden_size, bert_config.hidden_size))
            layers.append(torch.nn.ReLU())

        if get_embeddings:
            layers.pop()
        else:
            layers.append(torch.nn.Linear(bert_config.hidden_size, 1))

        self.head = torch.nn.Sequential(*layers)

    def forward(self, inputs, attention_mask):
        output = self.protbert(inputs, attention_mask=attention_mask)

        return self.head(output.pooler_output)


def create_model(config, get_embeddings=False):
    bert_config = BertConfig(
        hidden_size=config['network']['hidden_size'],
        num_hidden_layers=config['network']['hidden_layers'],
        num_attention_heads=config['network']['attn_heads'],
        hidden_dropout_prob=config['network']['dropout'],
        max_position_embeddings=config['network']['max_position_embeddings']
    )
    model = IDPBERT(bert_config, config, get_embeddings).to(config['device'])

    return model


def setup_training(config, model):
    criterion = torch.nn.MSELoss()
    optimizer = torch.optim.AdamW(model.parameters(), lr=config['optim']['lr'])
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode='max',
        factor=config['sch']['factor'],
        patience=config['sch']['patience']
    )

    return criterion, optimizer, scheduler
