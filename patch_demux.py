import re

with open('demux_logic.py', 'r') as f:
    content = f.read()

# 1. Update _adapter_only_config to allow single adapters
old_code = """    for layout in layout_rules.get("valid_layouts", []) or []:
        adapter_layout = [
            str(site)
            for site in layout
            if _site_role(str(site), layout_rules) == "adapter"
        ]
        if len(adapter_layout) >= 2 and adapter_layout not in adapter_layouts:
            adapter_layouts.append(adapter_layout)"""

new_code = """    for layout in layout_rules.get("valid_layouts", []) or []:
        adapter_layout = [
            str(site)
            for site in layout
            if _site_role(str(site), layout_rules) == "adapter"
        ]
        if len(adapter_layout) >= 1 and adapter_layout not in adapter_layouts:
            adapter_layouts.append(adapter_layout)"""

content = content.replace(old_code, new_code)

# 2. Update matcher_core to not require 2 adapters if only 1 is in layout
with open('matcher_core.py', 'r') as f:
    matcher_content = f.read()

matcher_content = matcher_content.replace(
    'if len(adapters) < 2:',
    'if len(adapters) < 1:'
)
matcher_content = matcher_content.replace(
    'Expected at least 2 adapter matches',
    'Expected at least 1 adapter match'
)

# 3. Update _build_allowed_adapter_pairs to allow singletons
old_pairs = """    for left, right in combinations(adapters, 2):
        first, second = sorted((left, right), key=_match_sort_key)
        if first.side is not None and first.side == second.side:
            continue
        if _adapter_pair_allowed(first, second, layout_rules):
            adapter_pairs.append((first, second))

    return adapter_pairs"""

# We'll represent a single adapter as (match, None) or just (match, match) but that's messy.
# Better to update the logic to handle singletons.

with open('demux_logic.py', 'w') as f:
    f.write(content)
with open('matcher_core.py', 'w') as f:
    f.write(matcher_content)
