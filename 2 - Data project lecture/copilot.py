import pandas as pd
import yfinance as yf
from IPython.display import display
import matplotlib.pyplot as plt
from datetime import datetime
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})


class StockData:
    def __init__(self, ticker):
        self.ticker = ticker
        self.data = None

        self.download_data()
        self.calculate_running_means()

        self.calculate_returns()

        if self.data is not None:
            display(self.data.describe())


    def download_data(self):
        self.data = yf.download(self.ticker)

    def calculate_running_means(self, window_sizes=[200, 50]):
        if self.data is None:
            print("No data available. Download data first.")
            return

        for window in window_sizes:
            self.data[f'Adj Close RM {window}'] = self.data['Adj Close'].rolling(window=window).mean()
        

        # Create a dummy variable for when the two moving averages cross
        self.data['Cross'] = 0
        cross_points = (self.data['Adj Close RM 200'] > self.data['Adj Close RM 50']).astype(int)
        self.data.loc[cross_points.diff() != 0, 'Cross'] = 1

    def calculate_returns(self):
        if self.data is None:
            print("No data available. Download data first.")
            return

        # Calculate cumulative returns
        self.data['Return'] = self.data['Adj Close'].pct_change()
        self.data['Cumulative Return'] = (1 + self.data['Return']).cumprod()

        # Calculate alternative cumulative return
        self.data['Position'] =(self.data['Adj Close RM 50'] > self.data['Adj Close RM 200']).astype(int)
        #self.data['Position'].fillna(method='ffill', inplace=True)
        self.data['Alternative Return'] = self.data['Return'] * self.data['Position']
        self.data['Alternative Cumulative Return'] = (1 + self.data['Alternative Return']).cumprod()


    def plot(self, start_year, end_year, variables=['Adj Close', 'Adj Close RM 200', 'Adj Close RM 50']):
        if self.data is None:
            print("No data available. Download data first.")
            return

        start_date = datetime(start_year, 1, 1)
        end_date = datetime(end_year, 12, 31)

        df = self.data.loc[start_date:end_date].copy()

        # Rebase the cumulative returns if they are to be plotted
        if 'Cumulative Return' in variables or 'Alternative Cumulative Return' in variables:
            df['Cumulative Return'] /= df['Cumulative Return'].iloc[0]
            df['Alternative Cumulative Return'] /= df['Alternative Cumulative Return'].iloc[0]

        fig, ax = plt.subplots(figsize=(12,8))
        df[variables].plot(ax=ax)

        # Add vertical lines at the points where the two moving averages cross
        cross_points = df[df['Cross'] == 1].index
        for cross_point in cross_points:
            ax.axvline(x=cross_point, color='r', linestyle='--')

        ax.set_title(f'{self.ticker} Plot')
        ax.set_xlabel('Date')
        ax.set_ylabel('Value')
        plt.show()



    def plot_data_(self, start_year, end_year):
        if self.data is None:
            print("No data available. Download data first.")
            return

        start_date = datetime(start_year, 1, 1)
        end_date = datetime(end_year, 12, 31)

        df = self.data.loc[start_date:end_date].copy()
        fig, ax = plt.subplots(figsize=(12,8))
        df[['Adj Close', 'Adj Close RM 200', 'Adj Close RM 50']].plot(ax=ax)
        
        # Add vertical lines at the points where the two moving averages cross
        cross_points = df[df['Cross'] == 1].index
        for cross_point in cross_points:
            ax.axvline(x=cross_point, color='r', linestyle='--')
            
        ax.set_title(f'{self.ticker} Adjusted Close Price and Running Means')
        ax.set_xlabel('Date')
        ax.set_ylabel('Price')
        plt.show()


    
    def plot_data_(self, start_year, end_year):
        if self.data is None:
            print("No data available. Download data first.")
            return

        start_date = datetime(start_year, 1, 1)
        end_date = datetime(end_year, 12, 31)

        df = self.data.loc[start_date:end_date].copy()
        fig, ax = plt.subplots(figsize=(12,8))
        df[['Adj Close', 'Adj Close RM 200', 'Adj Close RM 50']].plot(ax=ax)
        
        # Add vertical lines at the points where the two moving averages cross
        cross_points = df[df['Cross'] == 1].index
        for cross_point in cross_points:
            ax.axvline(x=cross_point, color='r', linestyle='--')
            
        ax.set_title(f'{self.ticker} Adjusted Close Price and Running Means')
        ax.set_xlabel('Date')
        ax.set_ylabel('Price')
        plt.show()


    def plot_returns(self, start_year, end_year):
        if self.data is None:
            print("No data available. Download data first.")
            return

        start_date = datetime(start_year, 1, 1)
        end_date = datetime(end_year, 12, 31)

        df = self.data.loc[start_date:end_date].copy()
        fig, ax = plt.subplots(figsize=(12,8))

        # Rebase the cumulative returns to the start year
        df['Cumulative Return'] /= df['Cumulative Return'].iloc[0]
        df['Alternative Cumulative Return'] /= df['Alternative Cumulative Return'].iloc[0]
        
        df[['Cumulative Return', 'Alternative Cumulative Return']].plot(ax=ax)

        # Add vertical lines at the points where the two moving averages cross
        cross_points = df[df['Cross'] == 1].index
        for cross_point in cross_points:
            ax.axvline(x=cross_point, color='r', linestyle='--')

        ax.set_title(f'{self.ticker} Cumulative Returns and Crossover Points')
        ax.set_xlabel('Date')
        ax.set_ylabel('Cumulative Return')
        plt.show()