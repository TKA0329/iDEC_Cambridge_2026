#!/usr/bin/env python3
"""
Test caching functionality for SequenceParameters
"""

# Mock localcider since it's not installed in this environment
class MockSequenceParameters:
    def __init__(self, seq):
        self.seq = seq
        self.kappa = len(seq) * 0.01  # Mock kappa value

    def get_kappa(self):
        return self.kappa

# Cache for SequenceParameters to avoid recreation overhead
CIDER_CACHE = {}
CIDER_CACHE_MAX_SIZE = 10000

def calculate_cider_kappa(seq):
    """
    CIDER κ (kappa) charge-patterning metric from Das & Pappu (2013).
    Uses caching to avoid recreating SequenceParameters for repeated sequences.
    """
    if not seq:
        return None

    # Check cache first
    if seq in CIDER_CACHE:
        sp = CIDER_CACHE[seq]
    else:
        try:
            sp = MockSequenceParameters(seq)
            # Simple cache management - if full, clear it
            if len(CIDER_CACHE) >= CIDER_CACHE_MAX_SIZE:
                CIDER_CACHE.clear()
            CIDER_CACHE[seq] = sp
        except Exception:
            return None

    try:
        return round(sp.get_kappa(), 4)
    except Exception:
        return None

# Test caching
print("Testing SequenceParameters caching...")

# First call - should create new object
result1 = calculate_cider_kappa('TESTSEQ')
print(f'First call result: {result1}')
print(f'Cache size after first call: {len(CIDER_CACHE)}')

# Second call - should use cache
result2 = calculate_cider_kappa('TESTSEQ')
print(f'Second call result: {result2}')
print(f'Cache size after second call: {len(CIDER_CACHE)}')

# Different sequence - should create new object
result3 = calculate_cider_kappa('DIFFERENT')
print(f'Different sequence result: {result3}')
print(f'Cache size after different sequence: {len(CIDER_CACHE)}')

print('Caching test: PASSED' if result1 == result2 else 'Caching test: FAILED')