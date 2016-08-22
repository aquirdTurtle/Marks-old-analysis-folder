import smtplib

from email.mime.text import MIMEText

email = MIMEText('This is an email!\nMark', "plain")

email['Subject'] = "Test"
email["From"] = "o.mark.brown@gmail.com"
email["To"] = "Mark.O.Brown@colorado.edu"

mail = smtplib.SMTP('smtp.gmail.com', 587)

mail.ehlo()

mail.starttls()

mail.login('o.mark.brown@gmail.com', 'ruFF#nar5=cUeL<')
print(email.as_string())
mail.sendmail(email["From"], email["To"], email.as_string())
mail.quit()
