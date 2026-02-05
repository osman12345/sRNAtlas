# Deploying sRNA-seq Analysis WebTool to Streamlit Cloud

This guide provides step-by-step instructions to deploy your sRNA-seq Analysis WebTool to Streamlit Community Cloud (free hosting).

---

## Prerequisites

Before you begin, make sure you have:
1. A **GitHub account** (free): https://github.com/signup
2. A **Streamlit Cloud account** (free): https://share.streamlit.io (sign up with GitHub)
3. **Git** installed on your computer

---

## Step 1: Create a GitHub Repository

### Option A: Using GitHub Web Interface

1. Go to https://github.com/new
2. Fill in the repository details:
   - **Repository name**: `sRNA-WebTool` (or any name you prefer)
   - **Description**: "Web-based small RNA-seq analysis tool"
   - **Visibility**: Choose **Public** (required for free Streamlit Cloud)
   - Do NOT initialize with README (we already have one)
3. Click **Create repository**
4. Keep this page open - you'll need the repository URL

### Option B: Using GitHub CLI

```bash
# Install GitHub CLI if not already installed
# macOS: brew install gh
# Linux: sudo apt install gh

# Authenticate with GitHub
gh auth login

# Create repository
gh repo create sRNA-WebTool --public --description "Web-based small RNA-seq analysis tool"
```

---

## Step 2: Prepare Your Local Repository

Open Terminal (macOS/Linux) or Command Prompt/PowerShell (Windows) and navigate to your project:

```bash
# Navigate to the project directory
cd /path/to/sRNA_WebTool

# Initialize git repository (if not already done)
git init

# Add all files to staging
git add .

# Create the first commit
git commit -m "Initial commit: sRNA-seq Analysis WebTool"
```

---

## Step 3: Push to GitHub

Replace `YOUR_USERNAME` with your GitHub username:

```bash
# Add the remote repository
git remote add origin https://github.com/YOUR_USERNAME/sRNA-WebTool.git

# Push to GitHub
git branch -M main
git push -u origin main
```

If prompted, enter your GitHub credentials. For authentication, you may need to:
- Use a **Personal Access Token** instead of password
- Or set up **SSH keys**

### Creating a Personal Access Token (if needed)

1. Go to https://github.com/settings/tokens
2. Click **Generate new token (classic)**
3. Give it a name like "sRNA WebTool"
4. Select scopes: `repo` (full control of repositories)
5. Click **Generate token**
6. Copy the token and use it as your password when pushing

---

## Step 4: Deploy to Streamlit Cloud

### 4.1 Sign in to Streamlit Cloud

1. Go to https://share.streamlit.io
2. Click **Sign in with GitHub**
3. Authorize Streamlit to access your GitHub account

### 4.2 Create New App

1. Click **New app** button
2. Fill in the deployment form:

   | Field | Value |
   |-------|-------|
   | **Repository** | `YOUR_USERNAME/sRNA-WebTool` |
   | **Branch** | `main` |
   | **Main file path** | `app/main.py` |

3. Click **Advanced settings** (optional but recommended):
   - **Python version**: 3.10 or 3.11
   - You can add secrets here if needed later

4. Click **Deploy!**

### 4.3 Wait for Deployment

- The deployment typically takes 3-5 minutes
- You'll see a log of the installation process
- Once complete, your app will be live at: `https://YOUR_APP_NAME.streamlit.app`

---

## Step 5: Verify Your Deployment

1. Open your app URL
2. Test each module:
   - ✅ Quality Control (upload test data)
   - ✅ Differential Expression (upload count matrix)
   - ✅ GO/Pathway Enrichment
   - ✅ Settings page
   - ✅ Help documentation

---

## Troubleshooting

### Common Issues and Solutions

#### Issue: "ModuleNotFoundError"
**Solution**: Make sure all dependencies are in `requirements.txt` and the file is in the root directory.

#### Issue: App crashes on large file upload
**Solution**: Streamlit Cloud has a 200MB upload limit. For larger files, consider:
- Compressing input files
- Processing data in chunks
- Using cloud storage integration

#### Issue: "No module named 'app'"
**Solution**: Ensure the main file path is correct: `app/main.py`

#### Issue: Build fails with memory error
**Solution**: Some packages may be too large. Check if you can use lighter alternatives.

#### Issue: pyDESeq2 installation fails
**Solution**: This is rare but can happen. Try:
1. Specify exact version: `pydeseq2==0.4.0`
2. Check Streamlit Cloud logs for specific errors

---

## Updating Your App

After deployment, any changes pushed to GitHub will automatically trigger a redeploy:

```bash
# Make changes to your code
# Then commit and push
git add .
git commit -m "Update: description of changes"
git push origin main
```

The app will automatically redeploy within a few minutes.

---

## Custom Domain (Optional)

Streamlit Cloud allows custom domains on paid plans. For the free tier, your app will be available at:
```
https://YOUR_APP_NAME.streamlit.app
```

---

## Resource Limits (Free Tier)

Streamlit Community Cloud free tier includes:
- **Memory**: 1GB RAM
- **Storage**: Limited (use for temporary files only)
- **Compute**: Shared resources
- **Uptime**: Apps sleep after inactivity, wake on access

For production use with heavy workloads, consider upgrading or using alternative hosting.

---

## Sharing Your App

Once deployed, share your app by:

1. **Direct link**: `https://YOUR_APP_NAME.streamlit.app`
2. **Badge in README**:
   ```markdown
   [![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://YOUR_APP_NAME.streamlit.app)
   ```
3. **QR Code**: Generate at https://www.qr-code-generator.com/

---

## Quick Reference Commands

```bash
# Clone your repo on another machine
git clone https://github.com/YOUR_USERNAME/sRNA-WebTool.git

# Check status
git status

# View remote
git remote -v

# Pull latest changes
git pull origin main

# Push changes
git push origin main
```

---

## Need Help?

- **Streamlit Documentation**: https://docs.streamlit.io
- **Streamlit Community**: https://discuss.streamlit.io
- **GitHub Issues**: Create an issue in your repository

---

## Summary Checklist

- [ ] Create GitHub account
- [ ] Create Streamlit Cloud account (via GitHub)
- [ ] Create new GitHub repository
- [ ] Initialize local git repo
- [ ] Push code to GitHub
- [ ] Deploy on Streamlit Cloud
- [ ] Test all features
- [ ] Share your app URL!

**Estimated time**: 15-30 minutes

---

*Your sRNA-seq Analysis WebTool is now live and accessible to researchers worldwide!*
