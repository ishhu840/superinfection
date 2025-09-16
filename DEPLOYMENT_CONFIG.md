# Deployment Configuration Guide

## Production Deployment Setup

### Environment Requirements
- Python 3.8+
- 8GB+ RAM recommended
- 2GB+ disk space for data
- Internet connection for data downloads

### Quick Start Commands

#### 1. Install Dependencies
```bash
pip install -r requirements.txt
```

#### 2. Run Main Dashboard
```bash
streamlit run real_data_dashboard.py --server.port 8504
```

#### 3. Alternative Dashboards
```bash
# Comprehensive dashboard
streamlit run src/dashboard/comprehensive_dashboard.py --server.port 8503

# Simple app
streamlit run app/main_simple.py --server.port 8502

# Introduction page
streamlit run introduction_page.py --server.port 8505
```

### Production Configuration

#### Streamlit Config (.streamlit/config.toml)
```toml
[server]
port = 8504
headless = true
enableCORS = false
enableXsrfProtection = false

[browser]
gatherUsageStats = false

[theme]
primaryColor = "#1f77b4"
backgroundColor = "#ffffff"
secondaryBackgroundColor = "#f0f2f6"
textColor = "#262730"
```

#### Docker Configuration
```dockerfile
FROM python:3.9-slim

WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

EXPOSE 8504
CMD ["streamlit", "run", "real_data_dashboard.py", "--server.port=8504", "--server.address=0.0.0.0"]
```

### File Structure (Optimized)
```
super-virus/
├── real_data_dashboard.py          # Main production dashboard
├── introduction_page.py            # Introduction and methodology
├── requirements.txt                # Dependencies
├── DEPLOYMENT_CONFIG.md            # This file
├── CLEANUP_RECOMMENDATIONS.md      # Cleanup documentation
├── app/                           # Alternative applications
│   ├── main.py                   # Full-featured app
│   └── main_simple.py            # Simplified app
├── src/                          # Core source code
├── examples/                     # Demo and utility scripts
├── data/                        # Data files
├── real_viral_results/          # Analysis results
└── docs/                       # Documentation
```

### Health Checks

#### Application Status
- Dashboard loads without errors: ✅
- All visualizations render: ✅
- Data loads successfully: ✅
- No import errors: ✅

#### Performance Optimization
- Removed unused files: ✅
- Organized file structure: ✅
- Validated all imports: ✅
- Cleaned up redundant code: ✅

### Monitoring

#### Key Metrics to Monitor
1. Dashboard load time
2. Memory usage
3. Data processing time
4. User interaction responsiveness

#### Logs Location
- Application logs: `logs/`
- Streamlit logs: `.streamlit/`

### Backup and Recovery

#### Critical Files to Backup
- `real_data_dashboard.py`
- `real_viral_results/`
- `data/`
- `requirements.txt`

### Security Considerations

#### Production Security
- Run behind reverse proxy (nginx/Apache)
- Enable HTTPS
- Implement rate limiting
- Regular security updates

#### Data Privacy
- No sensitive data in logs
- Secure data transmission
- Regular data cleanup

### Troubleshooting

#### Common Issues
1. **Import errors**: Check requirements.txt installation
2. **Data loading errors**: Verify file paths and permissions
3. **Visualization errors**: Check Plotly version compatibility
4. **Memory issues**: Increase system RAM or optimize data loading

#### Support
- Check logs in `logs/` directory
- Verify all dependencies in requirements.txt
- Test with sample data first

### Deployment Checklist

- [ ] Dependencies installed
- [ ] Data files accessible
- [ ] Dashboard loads without errors
- [ ] All visualizations working
- [ ] Performance acceptable
- [ ] Security measures in place
- [ ] Monitoring configured
- [ ] Backup strategy implemented

## Ready for Production ✅

The application has been cleaned up, optimized, and is ready for deployment with:
- Streamlined file structure
- Removed unused code
- Validated all functionality
- Comprehensive documentation