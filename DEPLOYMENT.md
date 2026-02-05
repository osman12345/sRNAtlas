# Deployment Guide for sRNA-seq Analysis WebTool

This guide covers multiple deployment options from simplest to most robust.

---

## Table of Contents
1. [Streamlit Community Cloud (Free)](#1-streamlit-community-cloud-free)
2. [Google Cloud Run](#2-google-cloud-run)
3. [DigitalOcean App Platform](#3-digitalocean-app-platform)
4. [AWS App Runner](#4-aws-app-runner)
5. [Heroku](#5-heroku)
6. [Self-Hosted Server](#6-self-hosted-server)
7. [Publishing as a Python Package](#7-publishing-as-a-python-package)

---

## 1. Streamlit Community Cloud (Free)

**Best for:** Demos, small labs, testing

### Steps:

1. **Create a GitHub repository**
   ```bash
   git init
   git add .
   git commit -m "Initial commit"
   git remote add origin https://github.com/YOUR_USERNAME/sRNA-WebTool.git
   git push -u origin main
   ```

2. **Go to [share.streamlit.io](https://share.streamlit.io)**

3. **Click "New app" and connect your GitHub repo**

4. **Configure:**
   - Repository: `YOUR_USERNAME/sRNA-WebTool`
   - Branch: `main`
   - Main file path: `app/main.py`

5. **Deploy!**

### Limitations:
- 1GB memory limit
- Apps sleep after inactivity
- No persistent storage
- R/Bowtie2 not available (Python-only analysis)

### For Streamlit Cloud, create this simplified version:

Create `.streamlit/config.toml`:
```toml
[server]
maxUploadSize = 500
enableXsrfProtection = true

[theme]
primaryColor = "#1E88E5"
backgroundColor = "#FFFFFF"
secondaryBackgroundColor = "#F0F2F6"
textColor = "#262730"
```

---

## 2. Google Cloud Run

**Best for:** Production use, auto-scaling, pay-per-use

### Prerequisites:
- Google Cloud account
- `gcloud` CLI installed

### Steps:

1. **Enable required APIs**
   ```bash
   gcloud services enable run.googleapis.com
   gcloud services enable containerregistry.googleapis.com
   ```

2. **Build and push Docker image**
   ```bash
   # Set your project ID
   PROJECT_ID="your-project-id"

   # Build image
   docker build -t gcr.io/$PROJECT_ID/srna-webtool .

   # Push to Google Container Registry
   docker push gcr.io/$PROJECT_ID/srna-webtool
   ```

3. **Deploy to Cloud Run**
   ```bash
   gcloud run deploy srna-webtool \
     --image gcr.io/$PROJECT_ID/srna-webtool \
     --platform managed \
     --region us-central1 \
     --memory 4Gi \
     --cpu 2 \
     --timeout 3600 \
     --allow-unauthenticated \
     --port 8501
   ```

4. **Get your URL**
   ```bash
   gcloud run services describe srna-webtool --format='value(status.url)'
   ```

### Cost Estimate:
- ~$0.00002400 per vCPU-second
- ~$0.00000250 per GiB-second
- Free tier: 2 million requests/month

---

## 3. DigitalOcean App Platform

**Best for:** Simple deployment, predictable pricing

### Steps:

1. **Push code to GitHub**

2. **Go to [DigitalOcean App Platform](https://cloud.digitalocean.com/apps)**

3. **Create App → Choose GitHub repo**

4. **Configure:**
   - Type: Web Service
   - Dockerfile path: `/Dockerfile`
   - HTTP Port: 8501
   - Instance Size: Basic ($12/mo) or Professional ($25/mo)

5. **Add environment variables (if needed)**

6. **Deploy**

### Create `app.yaml` for DigitalOcean:
```yaml
name: srna-webtool
services:
  - name: web
    dockerfile_path: Dockerfile
    http_port: 8501
    instance_count: 1
    instance_size_slug: professional-xs
    routes:
      - path: /
    health_check:
      http_path: /_stcore/health
```

---

## 4. AWS App Runner

**Best for:** AWS ecosystem users, auto-scaling

### Steps:

1. **Push image to ECR**
   ```bash
   # Create ECR repository
   aws ecr create-repository --repository-name srna-webtool

   # Get login
   aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin YOUR_ACCOUNT_ID.dkr.ecr.us-east-1.amazonaws.com

   # Tag and push
   docker tag srna-webtool:latest YOUR_ACCOUNT_ID.dkr.ecr.us-east-1.amazonaws.com/srna-webtool:latest
   docker push YOUR_ACCOUNT_ID.dkr.ecr.us-east-1.amazonaws.com/srna-webtool:latest
   ```

2. **Create App Runner service via Console or CLI**
   ```bash
   aws apprunner create-service \
     --service-name srna-webtool \
     --source-configuration '{
       "ImageRepository": {
         "ImageIdentifier": "YOUR_ACCOUNT_ID.dkr.ecr.us-east-1.amazonaws.com/srna-webtool:latest",
         "ImageRepositoryType": "ECR",
         "ImageConfiguration": {
           "Port": "8501"
         }
       }
     }' \
     --instance-configuration '{
       "Cpu": "2 vCPU",
       "Memory": "4 GB"
     }'
   ```

---

## 5. Heroku

**Best for:** Quick deployment, simple management

### Steps:

1. **Install Heroku CLI**
   ```bash
   # macOS
   brew tap heroku/brew && brew install heroku

   # Or download from heroku.com
   ```

2. **Create `heroku.yml`**
   ```yaml
   build:
     docker:
       web: Dockerfile
   run:
     web: streamlit run app/main.py --server.port=$PORT --server.address=0.0.0.0
   ```

3. **Deploy**
   ```bash
   heroku login
   heroku create srna-webtool
   heroku stack:set container
   git push heroku main
   ```

4. **Scale up (optional)**
   ```bash
   heroku ps:scale web=1:standard-2x  # $50/mo for more resources
   ```

### Limitations:
- Dynos sleep after 30 min inactivity (free/hobby tier)
- 512MB RAM on free tier

---

## 6. Self-Hosted Server

**Best for:** Institutions, full control, data privacy

### Option A: Ubuntu Server with Docker

```bash
# On your server
sudo apt-get update
sudo apt-get install -y docker.io docker-compose

# Clone your repo
git clone https://github.com/YOUR_USERNAME/sRNA-WebTool.git
cd sRNA-WebTool

# Run with Docker Compose
docker-compose up -d

# Set up Nginx reverse proxy (recommended)
sudo apt-get install nginx
```

Create `/etc/nginx/sites-available/srna-webtool`:
```nginx
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://localhost:8501;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 86400;
    }
}
```

Enable and get SSL:
```bash
sudo ln -s /etc/nginx/sites-available/srna-webtool /etc/nginx/sites-enabled/
sudo certbot --nginx -d your-domain.com
sudo systemctl restart nginx
```

### Option B: Using systemd (without Docker)

Create `/etc/systemd/system/srna-webtool.service`:
```ini
[Unit]
Description=sRNA-seq Analysis WebTool
After=network.target

[Service]
Type=simple
User=www-data
WorkingDirectory=/opt/sRNA-WebTool
Environment="PATH=/opt/sRNA-WebTool/venv/bin"
ExecStart=/opt/sRNA-WebTool/venv/bin/streamlit run app/main.py --server.port=8501 --server.headless=true
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

```bash
sudo systemctl daemon-reload
sudo systemctl enable srna-webtool
sudo systemctl start srna-webtool
```

---

## 7. Publishing as a Python Package

**Best for:** Distribution to bioinformatics community

### Create package structure:

```
srna_webtool/
├── pyproject.toml
├── README.md
├── LICENSE
├── src/
│   └── srna_webtool/
│       ├── __init__.py
│       ├── __main__.py
│       ├── app/
│       ├── modules/
│       ├── utils/
│       └── config/
```

### Create `pyproject.toml`:
```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "srna-webtool"
version = "1.0.0"
description = "Web-based small RNA-seq analysis tool"
readme = "README.md"
license = {text = "MIT"}
authors = [
    {name = "Your Name", email = "your.email@example.com"}
]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
]
requires-python = ">=3.9"
dependencies = [
    "streamlit>=1.28.0",
    "streamlit-option-menu>=0.3.6",
    "pandas>=2.0.0",
    "numpy>=1.24.0",
    "scipy>=1.10.0",
    "plotly>=5.15.0",
    "matplotlib>=3.7.0",
    "seaborn>=0.12.0",
    "biopython>=1.81",
    "scikit-learn>=1.3.0",
    "statsmodels>=0.14.0",
    "gprofiler-official",
    "pydeseq2",
]

[project.optional-dependencies]
full = ["pysam>=0.21.0", "rpy2>=3.5.0"]
dev = ["pytest", "black", "flake8"]

[project.scripts]
srna-webtool = "srna_webtool.__main__:main"

[project.urls]
Homepage = "https://github.com/YOUR_USERNAME/sRNA-WebTool"
Documentation = "https://github.com/YOUR_USERNAME/sRNA-WebTool#readme"
Repository = "https://github.com/YOUR_USERNAME/sRNA-WebTool"
```

### Create `src/srna_webtool/__main__.py`:
```python
import subprocess
import sys
from pathlib import Path

def main():
    app_path = Path(__file__).parent / "app" / "main.py"
    subprocess.run([sys.executable, "-m", "streamlit", "run", str(app_path)])

if __name__ == "__main__":
    main()
```

### Publish to PyPI:
```bash
# Build
pip install build twine
python -m build

# Upload to TestPyPI first
twine upload --repository testpypi dist/*

# Then to PyPI
twine upload dist/*
```

### Users can then install with:
```bash
pip install srna-webtool
srna-webtool  # Launches the app
```

---

## Comparison Table

| Method | Cost | Setup Time | Scalability | R/Bowtie2 Support |
|--------|------|------------|-------------|-------------------|
| Streamlit Cloud | Free | 5 min | Limited | ❌ |
| Google Cloud Run | ~$5-50/mo | 30 min | Excellent | ✅ |
| DigitalOcean | $12-25/mo | 20 min | Good | ✅ |
| AWS App Runner | ~$10-50/mo | 45 min | Excellent | ✅ |
| Heroku | $7-50/mo | 15 min | Good | ✅ |
| Self-hosted | Server cost | 1-2 hours | Full control | ✅ |
| PyPI Package | Free | 1 hour | N/A (local) | ✅ |

---

## Recommendations

1. **For demos/testing:** Streamlit Community Cloud
2. **For research labs:** Google Cloud Run or DigitalOcean
3. **For institutions:** Self-hosted with Docker
4. **For distribution:** PyPI package + Docker image

---

## Security Considerations

- Always use HTTPS in production
- Set up authentication if handling sensitive data
- Regularly update dependencies
- Consider data privacy regulations (GDPR, HIPAA)
- Use environment variables for secrets
