# Code Cleanup and Optimization Recommendations

## Files Identified for Cleanup

### Potentially Unused Files
1. **simple_test.py** - Basic test script that duplicates functionality in test_pipeline.py
2. **test_enhanced_pipeline.py** - May be redundant with test_pipeline.py
3. **analyze_downloaded_genomes.py** - Standalone script that may be integrated into main pipeline
4. **demo_data_collection.py** - Demo script that could be moved to examples/
5. **generate_enhanced_results.py** - Utility script that could be integrated

### Dashboard Files Analysis
- **real_data_dashboard.py** - Main production dashboard ✅ KEEP
- **src/dashboard/comprehensive_dashboard.py** - Alternative dashboard implementation
- **app/main.py** - Full-featured app ✅ KEEP
- **app/main_simple.py** - Simplified version for basic packages
- **introduction_page.py** - Standalone introduction page

### Recommendations

#### High Priority Cleanup
1. **Consolidate test files**: Merge simple_test.py functionality into test_pipeline.py
2. **Remove duplicate imports**: Clean up unused imports in real_data_dashboard.py
3. **Organize utility scripts**: Move demo and utility scripts to appropriate subdirectories

#### Medium Priority
1. **Dashboard consolidation**: Choose primary dashboard and archive alternatives
2. **Code deduplication**: Identify and merge similar functions across files
3. **Import optimization**: Remove unused imports throughout codebase

#### File Structure Optimization
```
super-virus/
├── src/                    # Core source code
├── app/                    # Streamlit applications
├── examples/               # Demo and example scripts
├── tests/                  # All test files
├── docs/                   # Documentation
├── data/                   # Data files
└── scripts/                # Utility scripts
```

## Import Analysis for real_data_dashboard.py

### Currently Used Imports ✅
- streamlit as st
- pandas as pd
- plotly.express as px
- plotly.graph_objects as go
- numpy as np (used in chord diagram)
- defaultdict, Counter (used in virus detection analysis)
- networkx as nx (used in network visualization)
- json, pathlib

### All imports are being used - no cleanup needed for imports

## Next Steps
1. Archive or remove unused test files
2. Consolidate dashboard implementations
3. Organize utility scripts into proper directories
4. Update documentation to reflect cleaned structure