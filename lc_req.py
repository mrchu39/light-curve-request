from flask import Flask, render_template, url_for, flash, redirect, request
from forms import ZTFName, UploadToFritz
import func, snid
import matplotlib.pyplot as plt
import sncosmo, os
import numpy as np

z = 0.060574239946858476
z_std = 0.023994157056121096
x1 = -0.14238796934437334
x1_std = 1.4557579021314682
c = 0.08928354223298558
c_std = 0.15670291093588692
x0 = 0.0007648532623426458
x0_std = 0.0004363803462578883

app = Flask(__name__)
app.config['SECRET_KEY'] = '55ce971fc21bc4d516c94d98bee2f5f8'

@app.route("/", methods=['GET', 'POST'])
def home():
    form = ZTFName()
    if form.validate_on_submit():
        resp = func.api('GET', func.BASEURL+'api/sources/'+form.ztfname.data)
        if resp['status'] == 'error':
            flash(resp['message'], 'danger')
            return redirect(url_for('home'))
        else:
            return redirect(url_for('query', ztfname=form.ztfname.data))
    return render_template('home.html', form=form)

@app.route("/how-to-use")
def how_to_use():
    return render_template('how_to_use.html')

@app.route("/query/<ztfname>", methods=['GET', 'POST'])
def query(ztfname):

    form = UploadToFritz()

    try:
        data, result, fitted_model = snid.model_lc(ztfname)
    except RuntimeError:
        flash('Runtime Error encountered, no good fit likely.', 'danger')
        return redirect(url_for('home'))

    x1_nstds = np.round(np.abs((result.parameters[3]-x1)/x1_std), 1)
    c_nstds = np.round(np.abs((result.parameters[4]-c))/c_std, 1)

    sncosmo.plot_lc(data, model=fitted_model)
    plt.savefig('static/temp.png')

    fit_info = {
        'ztfname' : ztfname,
        'x1_nstds' : x1_nstds,
        'c_nstds' : c_nstds,
        'peak_absmag' : np.round(snid.get_peak_absmag(result.parameters[0], result.parameters[2]), 1)
    }

    if form.validate_on_submit():
        resp = snid.post_lc(ztfname)

        if resp == False:
            flash(ztfname + ' LC is already up to date.', 'success')
            return redirect(url_for('query', ztfname=ztfname))
        if resp['status'] == 'error':
            flash(resp['message'], 'danger')
            return redirect(url_for('query', ztfname=ztfname))
        else:
            flash(f'{ztfname} LC successfully uploaded.', 'success')
            return redirect(url_for('query', ztfname=ztfname))

    return render_template('query.html', fit_info=fit_info, form=form)

if __name__ == '__main__':
    app.run(debug=True)
