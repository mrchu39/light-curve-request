from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired, Length

class ZTFName(FlaskForm):
    ztfname = StringField('ZTFName', validators=[DataRequired()])
    submit = SubmitField('Query Fritz')

class UploadToFritz(FlaskForm):
    submit = SubmitField('Upload to Fritz')
