FROM python:3.13-slim

# Set the working directory
WORKDIR /my-app

# Copy the requirements file
COPY requirements.txt .

# Install the dependencies
RUN pip install -r requirements.txt

# Copy the application code
COPY . .

# Specify the command to run
CMD ["python", "plot_interactive.py"]