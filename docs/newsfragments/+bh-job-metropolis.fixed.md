BasinHoppingJob::run applies exp(-dE/(kB T)) only when the hop is uphill and temperature is positive. Temperature at or below zero rejects that hop instead of dividing by temperature.
