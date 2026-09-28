In-process jobs poll a cancel token before each job and before each dispatch. ``cancel_state`` sets that token and returns 1. A compiled relax still finishes the current call.
