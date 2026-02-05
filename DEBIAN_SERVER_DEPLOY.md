# Deploying sRNA-seq Analysis WebTool on Debian Server

Complete guide for deploying on a remote Debian/Ubuntu server accessed via SSH.

---

## Prerequisites

- Debian 10+ or Ubuntu 20.04+ server
- SSH access with sudo privileges
- At least 2GB RAM (4GB+ recommended)
- Domain name (optional, for HTTPS)

---

## Step 1: Connect to Your Server

```bash
ssh username@your-server-ip
```

---

## Step 2: Update System and Install Dependencies

```bash
# Update package list
sudo apt update && sudo apt upgrade -y

# Install Python and essential tools
sudo apt install -y python3 python3-pip python3-venv git nginx supervisor

# Check Python version (should be 3.9+)
python3 --version
```

---

## Step 3: Create Application User (Recommended)

```bash
# Create dedicated user for the app
sudo useradd -m -s /bin/bash srna-app

# Switch to app user
sudo su - srna-app
```

---

## Step 4: Clone/Upload Your Application

### Option A: Clone from GitHub

```bash
cd ~
git clone https://github.com/YOUR_USERNAME/sRNA-WebTool.git
cd sRNA-WebTool
```

### Option B: Upload via SCP (from your local machine)

```bash
# Run this on your LOCAL machine
scp -r /path/to/sRNA_WebTool username@your-server-ip:~/
```

Then on the server:
```bash
cd ~/sRNA_WebTool
```

---

## Step 5: Set Up Python Environment

```bash
# Create virtual environment
python3 -m venv venv

# Activate it
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install dependencies
pip install -r requirements.txt
```

---

## Step 6: Test the Application

```bash
# Quick test (Ctrl+C to stop)
streamlit run app/main.py --server.port 8501 --server.headless true
```

If working, you should see:
```
You can now view your Streamlit app in your browser.
Network URL: http://YOUR_SERVER_IP:8501
```

---

## Step 7: Configure Systemd Service (Auto-start on Boot)

Exit back to your sudo user:
```bash
exit  # Exit srna-app user
```

Create systemd service file:

```bash
sudo nano /etc/systemd/system/srna-webtool.service
```

Paste this content:

```ini
[Unit]
Description=sRNA-seq Analysis WebTool
After=network.target

[Service]
Type=simple
User=srna-app
Group=srna-app
WorkingDirectory=/home/srna-app/sRNA-WebTool
Environment="PATH=/home/srna-app/sRNA-WebTool/venv/bin"
ExecStart=/home/srna-app/sRNA-WebTool/venv/bin/streamlit run app/main.py --server.port 8501 --server.headless true --server.address 0.0.0.0
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

Save and exit (Ctrl+X, Y, Enter).

Enable and start the service:

```bash
# Reload systemd
sudo systemctl daemon-reload

# Enable on boot
sudo systemctl enable srna-webtool

# Start the service
sudo systemctl start srna-webtool

# Check status
sudo systemctl status srna-webtool
```

---

## Step 8: Configure Nginx Reverse Proxy

This allows you to access the app on port 80/443 with proper WebSocket support.

```bash
sudo nano /etc/nginx/sites-available/srna-webtool
```

Paste this configuration:

```nginx
server {
    listen 80;
    server_name your-domain.com;  # Or use server IP

    # Increase max upload size (for FASTQ files)
    client_max_body_size 500M;

    location / {
        proxy_pass http://127.0.0.1:8501;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 86400;
    }

    # Streamlit specific WebSocket endpoints
    location /_stcore/stream {
        proxy_pass http://127.0.0.1:8501/_stcore/stream;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_read_timeout 86400;
    }
}
```

Enable the site:

```bash
# Create symlink
sudo ln -s /etc/nginx/sites-available/srna-webtool /etc/nginx/sites-enabled/

# Remove default site (optional)
sudo rm /etc/nginx/sites-enabled/default

# Test nginx config
sudo nginx -t

# Restart nginx
sudo systemctl restart nginx
```

---

## Step 9: Configure Firewall

```bash
# Allow HTTP and HTTPS
sudo ufw allow 80/tcp
sudo ufw allow 443/tcp

# If you want direct access to Streamlit port (optional)
sudo ufw allow 8501/tcp

# Enable firewall if not already
sudo ufw enable

# Check status
sudo ufw status
```

---

## Step 10: Add SSL/HTTPS (Recommended)

Using Let's Encrypt for free SSL:

```bash
# Install certbot
sudo apt install -y certbot python3-certbot-nginx

# Get certificate (replace with your domain)
sudo certbot --nginx -d your-domain.com

# Auto-renewal is set up automatically
# Test renewal
sudo certbot renew --dry-run
```

---

## Step 11: Access Your Application

Your app is now available at:

- **With domain**: `https://your-domain.com`
- **Without domain**: `http://your-server-ip`
- **Direct Streamlit port**: `http://your-server-ip:8501`

---

## Useful Commands

### Service Management

```bash
# Check status
sudo systemctl status srna-webtool

# Restart app
sudo systemctl restart srna-webtool

# Stop app
sudo systemctl stop srna-webtool

# View logs
sudo journalctl -u srna-webtool -f

# View last 100 lines
sudo journalctl -u srna-webtool -n 100
```

### Updating the Application

```bash
# Switch to app user
sudo su - srna-app
cd ~/sRNA-WebTool

# If using Git
git pull origin main

# Or upload new files via SCP, then:
source venv/bin/activate
pip install -r requirements.txt

# Exit and restart
exit
sudo systemctl restart srna-webtool
```

### Nginx Commands

```bash
# Test config
sudo nginx -t

# Restart
sudo systemctl restart nginx

# View logs
sudo tail -f /var/log/nginx/error.log
```

---

## Installing Full Pipeline Tools (Optional)

For the complete pipeline with alignment capabilities:

```bash
# Install system dependencies
sudo apt install -y build-essential zlib1g-dev libbz2-dev liblzma-dev

# Install Bowtie2
sudo apt install -y bowtie2

# Install Samtools
sudo apt install -y samtools

# Verify
bowtie2 --version
samtools --version
```

---

## Performance Tuning

### Increase File Upload Limit

Edit Streamlit config:

```bash
sudo su - srna-app
mkdir -p ~/.streamlit
nano ~/.streamlit/config.toml
```

Add:
```toml
[server]
maxUploadSize = 1000  # MB
maxMessageSize = 1000  # MB

[browser]
gatherUsageStats = false
```

### Increase System Limits

```bash
sudo nano /etc/security/limits.conf
```

Add:
```
srna-app soft nofile 65536
srna-app hard nofile 65536
```

---

## Troubleshooting

### App not starting

```bash
# Check logs
sudo journalctl -u srna-webtool -n 50

# Test manually
sudo su - srna-app
cd ~/sRNA-WebTool
source venv/bin/activate
streamlit run app/main.py --server.headless true
```

### WebSocket connection errors

Ensure nginx config has proper WebSocket headers (see Step 8).

### Permission errors

```bash
# Fix ownership
sudo chown -R srna-app:srna-app /home/srna-app/sRNA-WebTool
```

### Port already in use

```bash
# Find what's using the port
sudo lsof -i :8501

# Kill if needed
sudo kill -9 PID
```

### Out of memory

```bash
# Check memory
free -h

# Add swap if needed
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab
```

---

## Security Recommendations

1. **Keep system updated**: `sudo apt update && sudo apt upgrade`
2. **Use SSH keys** instead of passwords
3. **Enable fail2ban**: `sudo apt install fail2ban`
4. **Use HTTPS** with Let's Encrypt
5. **Regular backups** of user data
6. **Limit SSH access** to specific IPs if possible

---

## Quick Reference

| Task | Command |
|------|---------|
| Start app | `sudo systemctl start srna-webtool` |
| Stop app | `sudo systemctl stop srna-webtool` |
| Restart app | `sudo systemctl restart srna-webtool` |
| View logs | `sudo journalctl -u srna-webtool -f` |
| Update app | `git pull && sudo systemctl restart srna-webtool` |
| Check nginx | `sudo nginx -t && sudo systemctl restart nginx` |

---

## Summary

Your sRNA-seq Analysis WebTool is now running on your Debian server with:

- ✅ Automatic startup on boot
- ✅ Nginx reverse proxy
- ✅ SSL/HTTPS support (if configured)
- ✅ Large file upload support
- ✅ Production-ready configuration

Access it at: `http://your-server-ip` or `https://your-domain.com`
